// mathilda_interpreter.cpp — see mathilda_interpreter.hpp.
//
// Maps Mathilda's NDJSON cell events (mathilda_ffi_eval_cell) onto Jupyter
// display messages:
//   stream  -> publish_stream("stdout")           (Print output)
//   message -> publish_stream("stderr")           (Head::tag warnings)
//   expr    -> {text/latex: "$..$", text/plain}   (the result; see below)
//   plot    -> {application/vnd.plotly.v1+json}   (Graphics / Plot)
//   image   -> {application/json}                 (Image; PNG is a follow-up)
//   usage   -> {text/plain}                       (?name docstring)
//   names   -> {text/plain}                        (?Pat* match list)
//   error   -> publish_execution_error            (parse error)
//
// A cell may hold several statements, so several "expr" events can arrive.
// Mathematica shows every non-suppressed result; Jupyter reserves a single
// execute_result (Out[n]) per cell. We publish every intermediate result as
// display_data and the LAST one as execute_result, preserving interleaving
// with Print/message output by flushing the pending result just before the
// next event.

#include "mathilda_interpreter.hpp"

#include <functional>
#include <string>
#include <vector>

#include <xeus/xhelper.hpp>

extern "C" {
#include "mathilda_ffi.h"
}

namespace mathilda_kernel
{
    namespace
    {
        // Trampoline: the C ABI takes a plain function pointer + void* ctx; we
        // route it to a capturing std::function so the handler can call this
        // interpreter's (protected) publish_* members.
        using line_handler = std::function<void(const char*)>;

        void sink_trampoline(void* ctx, const char* json_line)
        {
            (*static_cast<line_handler*>(ctx))(json_line);
        }

        // Collect every event line (used by inspect, which runs a `?name` cell).
        void collect_trampoline(void* ctx, const char* json_line)
        {
            auto* out = static_cast<std::vector<std::string>*>(ctx);
            out->emplace_back(json_line);
        }

        bool is_symbol_char(char c)
        {
            return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')
                || (c >= '0' && c <= '9') || c == '$';
        }
    }

    void interpreter::configure_impl()
    {
        // Idempotent; honours $MATHILDA_HOME (set in kernel.json) for the
        // internal/ module tree, else falls back to the resolver ladder.
        mathilda_ffi_init();
    }

    void interpreter::execute_request_impl(send_reply_callback cb,
                                           int execution_counter,
                                           const std::string& code,
                                           xeus::execute_request_config config,
                                           nl::json /*user_expressions*/)
    {
        const bool silent = config.silent;

        bool has_pending = false;
        nl::json pending;          // the mime bundle of the not-yet-flushed result
        bool had_error = false;
        std::string last_error;

        auto flush = [&](bool as_result) {
            if (!has_pending) return;
            has_pending = false;
            if (silent) return;
            if (as_result)
                publish_execution_result(execution_counter, pending, nl::json::object());
            else
                display_data(pending, nl::json::object(), nl::json::object());
        };

        line_handler handler = [&](const char* raw) {
            nl::json ev;
            try { ev = nl::json::parse(raw); }
            catch (...) { return; }
            const std::string type = ev.value("type", std::string());

            if (type == "done") { flush(true); return; }

            if (type == "expr") {
                flush(false);                       // previous result -> display_data
                nl::json data;
                data["text/plain"] = ev.value("payload", std::string());
                const std::string tex = ev.value("latex", std::string());
                if (!tex.empty()) data["text/latex"] = "$" + tex + "$";
                pending = std::move(data);
                has_pending = true;
                return;
            }

            flush(false);                           // keep output in source order

            if (type == "stream") {
                if (!silent) publish_stream("stdout", ev.value("text", std::string()));
            } else if (type == "message") {
                if (!silent) publish_stream("stderr", ev.value("text", std::string()) + "\n");
            } else if (type == "usage") {
                if (!silent) {
                    nl::json d; d["text/plain"] = ev.value("payload", std::string());
                    display_data(d, nl::json::object(), nl::json::object());
                }
            } else if (type == "names") {
                if (!silent && ev.contains("payload") && ev["payload"].is_array()) {
                    std::string joined;
                    for (const auto& n : ev["payload"]) {
                        if (!joined.empty()) joined += "   ";
                        joined += n.get<std::string>();
                    }
                    nl::json d; d["text/plain"] = joined;
                    display_data(d, nl::json::object(), nl::json::object());
                }
            } else if (type == "plot") {
                if (!silent && ev.contains("payload")) {
                    nl::json d;
                    d["application/vnd.plotly.v1+json"] = ev["payload"];
                    d["text/plain"] = "-Graphics-";
                    display_data(d, nl::json::object(), nl::json::object());
                }
            } else if (type == "image") {
                if (!silent && ev.contains("payload")) {
                    nl::json d;
                    d["application/json"] = ev["payload"];
                    d["text/plain"] = "-Image-";
                    display_data(d, nl::json::object(), nl::json::object());
                }
            } else if (type == "error") {
                had_error = true;
                last_error = ev.value("message", std::string());
                if (!silent)
                    publish_execution_error("Error", last_error,
                                            std::vector<std::string>{ last_error });
            }
            // "line": the kernel's $Line; Jupyter keeps its own count — ignore.
        };

        mathilda_ffi_eval_cell(code.c_str(), execution_counter, /*cell=*/1,
                               sink_trampoline, &handler);

        if (had_error)
            cb(xeus::create_error_reply("Error", last_error));
        else
            cb(xeus::create_successful_reply());
    }

    nl::json interpreter::complete_request_impl(const std::string& code, int cursor_pos)
    {
        // The identifier immediately left of the cursor is the prefix to match.
        int start = cursor_pos;
        while (start > 0 && is_symbol_char(code[start - 1])) --start;
        const std::string prefix = code.substr(start, cursor_pos - start);

        nl::json matches = nl::json::array();
        char* raw = mathilda_ffi_complete(prefix.c_str());
        if (raw) {
            try { matches = nl::json::parse(raw); }
            catch (...) { matches = nl::json::array(); }
            mathilda_ffi_free(raw);
        }
        return xeus::create_complete_reply(matches, start, cursor_pos);
    }

    nl::json interpreter::inspect_request_impl(const std::string& code,
                                               int cursor_pos,
                                               int /*detail_level*/)
    {
        // Widen to the whole identifier under the cursor, then ask `?name`.
        int start = cursor_pos, end = cursor_pos;
        while (start > 0 && is_symbol_char(code[start - 1])) --start;
        while (end < static_cast<int>(code.size()) && is_symbol_char(code[end])) ++end;
        const std::string sym = code.substr(start, end - start);
        if (sym.empty()) return xeus::create_inspect_reply(false);

        std::vector<std::string> events;
        const std::string query = "?" + sym;
        mathilda_ffi_eval_cell(query.c_str(), 0, /*cell=*/1,
                               collect_trampoline, &events);

        for (const auto& line : events) {
            nl::json ev;
            try { ev = nl::json::parse(line); } catch (...) { continue; }
            if (ev.value("type", std::string()) == "usage") {
                nl::json data; data["text/plain"] = ev.value("payload", std::string());
                return xeus::create_inspect_reply(true, data);
            }
        }
        return xeus::create_inspect_reply(false);
    }

    nl::json interpreter::is_complete_request_impl(const std::string& code)
    {
        return xeus::create_is_complete_reply(
            mathilda_ffi_is_complete(code.c_str()) ? "complete" : "incomplete");
    }

    nl::json interpreter::kernel_info_request_impl()
    {
        const std::string version = mathilda_ffi_version();
        const std::string banner =
            "Mathilda " + version + " — a Mathematica-like computer algebra system.";
        return xeus::create_info_reply(
            /*implementation*/        "xmathilda",
            /*implementation_version*/version,
            /*language_name*/         "mathilda",
            /*language_version*/      version,
            /*language_mimetype*/     "text/x-mathematica",
            /*language_file_extension*/".m",
            /*pygments_lexer*/        "mathematica",
            /*codemirror_mode*/       std::string("mathematica"),
            /*nbconvert_exporter*/    "",
            /*banner*/                banner);
    }

    nl::json interpreter::shutdown_request_impl(bool /*restart*/)
    {
        // The evaluator's state is process-global; a restart is a fresh process,
        // so there is nothing to tear down here.
        return nl::json::object();
    }

    nl::json interpreter::interrupt_request_impl()
    {
        // Mathilda's evaluator is not yet interruptible mid-computation (an abort
        // flag checked in the eval loop is a follow-up); acknowledge the request.
        return nl::json::object();
    }
}
