# Mac App Store distribution runbook

This document describes how to build, sign, and upload the Mathilda Notebook
(the Tauri v2 app in `frontend/`) to the Mac App Store (MAS).

Reference: <https://v2.tauri.app/distribute/app-store/>

> **Scope / status**
>
> - **arm64-only** for now. The bundled math libraries (GMP, MPFR, pcre2) are
>   copied from the local arm64 Homebrew install; building x86_64 versions of
>   these here is not feasible. Universal (arm64 + x86_64) is **future work** —
>   it requires x86_64 builds of GMP/MPFR/pcre2/raylib and `lipo`-merged
>   dylibs + a universal engine binary.
> - All four Homebrew dylibs are bundled (gmp, mpfr, pcre2, raylib). raylib is
>   kept even though the notebook renders graphics via Plotly JSON, to avoid any
>   regression risk in the engine's `USE_GRAPHICS` code paths.

---

## 0. What is automated vs. what needs YOUR Apple account

| Step | Automatable in-repo | Requires your Apple Developer account |
|------|:-------------------:|:-------------------------------------:|
| Self-contained sidecar (`build-sidecar-appstore.sh`) | ✅ | — |
| Entitlements / app-store Tauri config | ✅ | — |
| Apple Developer Program enrollment | — | ✅ |
| Register App ID + App Store Connect record | — | ✅ |
| Create Distribution + Installer certificates | — | ✅ |
| Create + download provisioning profile | — | ✅ |
| Signed `.app` / `.pkg` build | ✅ (once certs+profile are installed) | needs your identity |
| Upload to App Store Connect | — | ✅ (API key / Apple ID) |

Everything marked "requires your Apple Developer account" **cannot be done in
this repo** — it needs a paid Apple Developer account, private signing keys, and
web actions on developer.apple.com / App Store Connect.

---

## 1. Self-contained sidecar (already solved in-repo)

MAS apps run inside the **App Sandbox**, which cannot load dylibs from
`/opt/homebrew`. The engine is spawned as a Tauri **sidecar**
(`frontend/src-tauri/binaries/mathilda-aarch64-apple-darwin`) and normally links
Homebrew dylibs. `frontend/build-sidecar-appstore.sh` fixes this:

```sh
cd frontend
./build-sidecar-appstore.sh
```

It builds the engine, copies the sidecar, bundles the required Homebrew dylibs
(and any transitive Homebrew deps) into `src-tauri/binaries/libs/`, rewrites all
install names to `@rpath`/`@loader_path`, adds `LC_RPATH`s
(`@executable_path/../Frameworks` and `@executable_path/libs`), and ad-hoc
re-signs everything. It verifies at the end that **no `/opt/homebrew` paths
remain**.

`tauri.appstore.conf.json` lists those dylibs under `bundle.macOS.frameworks`,
so Tauri copies them into `Mathilda.app/Contents/Frameworks/` at bundle time and
includes them in the final signature. The sidecar (in `Contents/MacOS/`)
resolves them via `@rpath` → `../Frameworks`.

---

## 2. Apple Developer enrollment (requires your account)

1. Enroll in the **Apple Developer Program** ($99/yr):
   <https://developer.apple.com/programs/enroll/>.
2. Note your **Team ID** (10 chars, e.g. `A1B2C3D4E5`) — Membership page.

---

## 3. Register the App ID and App Store Connect record (requires your account)

1. developer.apple.com → **Certificates, Identifiers & Profiles → Identifiers →
   +** → **App IDs → App**.
2. **Bundle ID** must EXACTLY equal `tauri.conf.json > identifier`, currently
   **`com.mathilda.notebook`**. (This is already a real reverse-DNS id — no fix
   needed.)
3. App Store Connect → **Apps → +** → **New App**, select macOS, and pick the
   same bundle id.

---

## 4. Certificates (requires your account)

Create and install (double-click the downloaded `.cer` to add to Keychain):

- **Apple Distribution** — signs the `.app`.
  Identity string: `Apple Distribution: YOUR NAME (TEAMID)`
- **Mac Installer Distribution** — signs the `.pkg` for upload.
  Identity string: `3rd Party Mac Developer Installer: YOUR NAME (TEAMID)`

Verify they are present:

```sh
security find-identity -v -p codesigning     # should list "Apple Distribution: ..."
security find-identity -v | grep "Installer"  # should list "3rd Party Mac Developer Installer: ..."
```

---

## 5. Provisioning profile (requires your account)

1. developer.apple.com → **Profiles → +** → **Mac App Store** distribution.
2. Select the `com.mathilda.notebook` App ID and the Apple Distribution cert.
3. Download the profile and copy it to the exact filename referenced by
   `tauri.appstore.conf.json`:

   ```sh
   cp ~/Downloads/Mathilda_MAS.provisionprofile \
      frontend/src-tauri/embedded.provisionprofile
   ```

   Tauri copies this to `Mathilda.app/Contents/embedded.provisionprofile`, where
   macOS expects it.

---

## 6. Fill in the entitlements

`frontend/src-tauri/Entitlements.plist` ships with placeholders. Replace them:

- `$IDENTIFIER` → `com.mathilda.notebook` (must match `tauri.conf.json`)
- `$TEAM_ID` → your 10-char Team ID

Entitlements included (and why):

| Entitlement | Value | Why |
|-------------|-------|-----|
| `com.apple.security.app-sandbox` | `true` | **Mandatory** for MAS. |
| `com.apple.application-identifier` | `com.mathilda.notebook` | Must match bundle id + profile. |
| `com.apple.developer.team-identifier` | `<TEAMID>` | Must match cert + profile. |
| `com.apple.security.files.user-selected.read-write` | `true` | The notebook opens/saves `.lb` files via the dialog plugin; the sandbox only allows user-picked files with this entitlement. |

**Not included:** `com.apple.security.network.*` — the engine talks to the app
over local stdio pipes only; there is no network use. Add `network.client` only
if a real network feature is introduced.

---

## 7. Build the signed `.app`

```sh
cd frontend
./build-sidecar-appstore.sh          # (re)build the self-contained sidecar

export APPLE_SIGNING_IDENTITY="Apple Distribution: YOUR NAME (TEAMID)"
npm run tauri build -- --bundles app --config src-tauri/tauri.appstore.conf.json
```

This produces `frontend/src-tauri/target/release/bundle/macos/Mathilda.app`,
signed with the App Sandbox entitlements and the embedded provisioning profile.

> The `--config src-tauri/tauri.appstore.conf.json` overlay is merged onto the
> base `tauri.conf.json`. Normal `.dmg` builds (`npm run tauri build`) do **not**
> use it and are unaffected.

---

## 8. Wrap into a signed `.pkg`

MAS uploads require a `.pkg` signed with the **installer** identity:

```sh
xcrun productbuild \
  --sign "3rd Party Mac Developer Installer: YOUR NAME (TEAMID)" \
  --component "src-tauri/target/release/bundle/macos/Mathilda.app" /Applications \
  Mathilda.pkg
```

---

## 9. Upload to App Store Connect (requires your account)

Using an App Store Connect API key (recommended):

```sh
xcrun altool --upload-app --type macos --file Mathilda.pkg \
  --apiKey <KEY_ID> --apiIssuer <ISSUER_ID>
```

(The newer `xcrun notarytool` is for notarization of Developer-ID/`.dmg`
distribution, not MAS uploads; MAS apps are reviewed, not notarized.)

Then finish submission in App Store Connect (screenshots, metadata, review
notes) and submit for review.

---

## 10. Common pitfalls

- **"App sandbox not enabled" rejection.** Entitlements are applied at
  *signing* time, not build time. You must build with signing enabled
  (`APPLE_SIGNING_IDENTITY` set) and the `tauri.appstore.conf.json` overlay so
  `bundle.macOS.entitlements` points at `Entitlements.plist`. A plain unsigned
  build will be rejected. Verify:

  ```sh
  codesign -d --entitlements - \
    src-tauri/target/release/bundle/macos/Mathilda.app
  ```

  should show `com.apple.security.app-sandbox = true`.

- **Root-only-readable files (installer error 409 / "package ... not
  readable").** Every file in the bundle must be world-readable and not
  owned by root. `build-sidecar-appstore.sh` `chmod u+w`s the copied dylibs;
  if you copy files manually, ensure `chmod -R a+r` and non-root ownership.

- **Architecture mismatch.** This pipeline is arm64-only. Uploading an
  arm64-only build is accepted, but the app will not run on Intel Macs. For
  universal support, build x86_64 GMP/MPFR/pcre2/raylib, `lipo`-merge each dylib
  and the engine binary, and re-run the install-name rewrites.

- **`/opt/homebrew` sneaking back in.** After any engine rebuild, re-run
  `build-sidecar-appstore.sh` — it re-verifies no Homebrew paths remain. Confirm
  manually with:

  ```sh
  otool -L src-tauri/binaries/mathilda-aarch64-apple-darwin | grep homebrew || echo clean
  ```

- **Signature invalidated after `install_name_tool`.** Any post-signing edit to
  a Mach-O invalidates its signature. The script re-signs after rewriting; the
  real Tauri build re-signs everything again with your Distribution identity, so
  order matters: build sidecar → tauri build (signs) → productbuild.
