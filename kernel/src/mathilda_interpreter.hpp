// mathilda_interpreter.hpp — xeus Jupyter kernel for the Mathilda CAS.
//
// A thin C++ adaptor: xeus implements the Jupyter wire protocol; this class
// forwards each request to Mathilda's in-process C ABI (src/ffi/mathilda_ffi.h)
// and maps the kernel's NDJSON events onto Jupyter display messages. No Python
// is involved at runtime — only the C evaluator and xeus.
//
// One kernel == one process: the Mathilda evaluator is a process-global
// singleton (its symbol table is global and not reentrant), which is exactly
// Jupyter's one-kernel-per-process model.
#ifndef MATHILDA_KERNEL_INTERPRETER_HPP
#define MATHILDA_KERNEL_INTERPRETER_HPP

#include <string>

#include <nlohmann/json.hpp>
#include <xeus/xeus.hpp>            /* XEUS_VERSION_MAJOR */
#include <xeus/xinterpreter.hpp>

// xeus 6 reshaped two of xinterpreter's pure virtuals, so the kernel has to pick
// a branch at compile time (issue #85 — Sage builds against xeus 5.2.x):
//
//              xeus 5.0 – 5.2.x                 xeus 6.0+
//   shutdown   void shutdown_request_impl()     nl::json shutdown_request_impl(bool restart)
//   interrupt  -- absent, the core never asks   nl::json interrupt_request_impl()
//
// A third difference is silent rather than fatal and is handled in the .cpp:
// create_info_reply() lost its leading `protocol_version` parameter in xeus 6,
// so the same all-std::string call binds under xeus 5 with every argument
// shifted one slot. Everything else the kernel uses — execute_request_impl and
// its config, complete/inspect/is_complete, the publish_* family, and all of
// main.cpp — is source-identical across 5.x and 6.x.
//
// This macro is the ONLY version test; add new divergences here, not inline.
// xeus < 5 is not supported (its execute_request_impl took an xrequest_context);
// CMakeLists.txt enforces that floor at configure time.
#if defined(XEUS_VERSION_MAJOR) && XEUS_VERSION_MAJOR >= 6
#  define MATHILDA_XEUS_6 1
#else
#  define MATHILDA_XEUS_6 0
#endif

namespace mathilda_kernel
{
    namespace nl = nlohmann;

    class interpreter : public xeus::xinterpreter
    {
    public:

        interpreter() = default;
        virtual ~interpreter() = default;

    private:

        void configure_impl() override;

        void execute_request_impl(send_reply_callback cb,
                                  int execution_counter,
                                  const std::string& code,
                                  xeus::execute_request_config config,
                                  nl::json user_expressions) override;

        nl::json complete_request_impl(const std::string& code,
                                       int cursor_pos) override;

        nl::json inspect_request_impl(const std::string& code,
                                      int cursor_pos,
                                      int detail_level) override;

        nl::json is_complete_request_impl(const std::string& code) override;

        nl::json kernel_info_request_impl() override;

#if MATHILDA_XEUS_6
        nl::json shutdown_request_impl(bool restart) override;

        nl::json interrupt_request_impl() override;
#else
        void shutdown_request_impl() override;   /* xeus 5: no flag, no reply */
#endif
    };
}

#endif
