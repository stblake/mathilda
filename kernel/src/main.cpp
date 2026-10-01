// main.cpp — entry point for the Mathilda xeus kernel (xmathilda).
//
// Boilerplate from the xeus kernel-authoring guide: load the connection file
// Jupyter passes as `-f <file>`, create a ZMQ context and the interpreter,
// then run the kernel. xeus owns the Jupyter wire protocol; all CAS behaviour
// is in mathilda_kernel::interpreter.

#include <memory>
#include <string>

#include <xeus/xkernel.hpp>
#include <xeus/xkernel_configuration.hpp>
#include <xeus/xhelper.hpp>

#include <xeus-zmq/xzmq_context.hpp>
#include <xeus-zmq/xserver_zmq.hpp>

#include "mathilda_interpreter.hpp"

int main(int argc, char* argv[])
{
    // `-f <connection_file>` — Jupyter's launch convention (removed from argv).
    const std::string connection_filename = xeus::extract_filename(argc, argv);

    auto context = xeus::make_zmq_context();
    auto interpreter = std::make_unique<mathilda_kernel::interpreter>();

    xeus::xconfiguration config = xeus::load_configuration(connection_filename);

    xeus::xkernel kernel(config,
                         xeus::get_user_name(),
                         std::move(context),
                         std::move(interpreter),
                         xeus::make_xserver_default);
    kernel.start();
    return 0;
}
