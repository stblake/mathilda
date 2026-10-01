#!/bin/sh
# MathildaXeus.command — double-click to open Mathilda in a Jupyter (xeus) front end.
#
# It lives at the REPOSITORY ROOT, beside MathildaNotebook.command, and resolves
# everything from its own location, so it works however it is invoked and
# wherever the repo is checked out.
#
# WHY THE .command SUFFIX. It is what makes the file clickable: Finder runs a
# .command file in Terminal (a plain executable opens in a text editor instead),
# and hides the extension, so it reads as "MathildaXeus" in the folder.
#
# WHAT IT LAUNCHES. The xeus-based Jupyter kernel (kernel/, built as `xmathilda`)
# running under JupyterLab. Unlike the Tauri desktop notebook (MathildaNotebook),
# this is the standard Jupyter ecosystem front end — useful for users who already
# live in JupyterLab / VS Code / Colab. See kernel/README.md.
#
# WHERE THE TOOLCHAIN COMES FROM. The kernel needs the xeus toolchain + Jupyter,
# which are not vendored. Point this at a conda-forge env (recommended) via any
# of, in order: $MATHILDA_XEUS_ENV, an active $CONDA_PREFIX, or a conda/mamba/
# micromamba env named `mathilda-xeus`. Create one once with:
#
#   micromamba create -n mathilda-xeus -c conda-forge \
#       xeus xeus-zmq cppzmq nlohmann_json xtl cmake jupyterlab jupyter_client
#
# This script then builds + installs the kernel into that env if needed and
# opens JupyterLab with the Mathilda kernel available.

set -e

ROOT=$(cd "$(dirname "$0")" && pwd)

if [ "$(uname -s)" != "Darwin" ]; then
    echo "This launcher is macOS-only (it opens a Terminal via Finder)."
    echo "On Linux, activate your xeus env and run:"
    echo "    cmake -S '$ROOT/kernel' -B '$ROOT/kernel/build' -DCMAKE_PREFIX_PATH=\"\$CONDA_PREFIX\""
    echo "    cmake --build '$ROOT/kernel/build' -j && cmake --install '$ROOT/kernel/build' --prefix \"\$CONDA_PREFIX\""
    echo "    jupyter lab"
    exit 1
fi

PATH="$HOME/.cargo/bin:/opt/homebrew/bin:/usr/local/bin:/opt/local/bin:$PATH"
export PATH

# --- Locate the xeus / Jupyter environment ---------------------------------
ENV_PREFIX=""
if [ -n "$MATHILDA_XEUS_ENV" ] && [ -x "$MATHILDA_XEUS_ENV/bin/jupyter" ]; then
    ENV_PREFIX="$MATHILDA_XEUS_ENV"
elif [ -n "$CONDA_PREFIX" ] && [ -x "$CONDA_PREFIX/bin/jupyter" ]; then
    ENV_PREFIX="$CONDA_PREFIX"
else
    for mgr in micromamba mamba conda; do
        if command -v "$mgr" >/dev/null 2>&1; then
            # `<mgr> run -n mathilda-xeus` resolves the env without needing its
            # shell hook to be initialised in this non-interactive Terminal.
            CAND=$("$mgr" run -n mathilda-xeus printenv CONDA_PREFIX 2>/dev/null || true)
            if [ -n "$CAND" ] && [ -x "$CAND/bin/jupyter" ]; then
                ENV_PREFIX="$CAND"
                break
            fi
        fi
    done
fi

if [ -z "$ENV_PREFIX" ]; then
    echo "No xeus/Jupyter environment found."
    echo
    echo "Create one once (conda-forge has the whole toolchain, no Python kernel"
    echo "dependency at run time):"
    echo
    echo "    micromamba create -n mathilda-xeus -c conda-forge \\"
    echo "        xeus xeus-zmq cppzmq nlohmann_json xtl cmake jupyterlab jupyter_client"
    echo
    echo "then double-click this file again (or set MATHILDA_XEUS_ENV to its path)."
    exit 1
fi

JUPYTER="$ENV_PREFIX/bin/jupyter"
CMAKE="$ENV_PREFIX/bin/cmake"
[ -x "$CMAKE" ] || CMAKE=cmake
echo "Using environment: $ENV_PREFIX"

# --- Build + install the kernel if its kernelspec is missing ---------------
# The kernelspec, not the binary, is the source of truth: `jupyter` launches the
# kernel through it. (Re)build when absent; the Mathilda makefile builds
# libmathilda.a with its own GCC toolchain, independent of this env's compiler.
if ! JUPYTER_PATH="$ENV_PREFIX/share/jupyter" "$JUPYTER" kernelspec list 2>/dev/null | grep -q "xmathilda"; then
    echo "The Mathilda (xmathilda) kernel is not installed in this environment."
    echo "Building and installing it now — the first build takes a few minutes."
    echo
    "$CMAKE" -S "$ROOT/kernel" -B "$ROOT/kernel/build" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_PREFIX_PATH="$ENV_PREFIX" \
        -DCMAKE_INSTALL_PREFIX="$ENV_PREFIX" \
        -DCMAKE_INSTALL_RPATH="$ENV_PREFIX/lib"
    "$CMAKE" --build "$ROOT/kernel/build" -j
    "$CMAKE" --install "$ROOT/kernel/build"
    echo
fi

echo "Opening JupyterLab — pick the \"Mathilda\" kernel for a new notebook."
exec "$JUPYTER" lab
