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
#include <xeus/xinterpreter.hpp>

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

        nl::json shutdown_request_impl(bool restart) override;

        nl::json interrupt_request_impl() override;
    };
}

#endif
