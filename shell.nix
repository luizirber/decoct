let
  sources = import ./nix/sources.nix;
  rustPlatform = import ./nix/rust.nix { inherit sources; };
  pkgs = import sources.nixpkgs { overlays = [ (import sources.rust-overlay) ]; };
  mach-nix = import sources.mach-nix {
    pkgs = pkgs;
    python = "python39";
  };

  customPython = mach-nix.mkPython {
    requirements = ''
      sourmash>=4
      setuptools
    '';
  };
in
  with pkgs;

  pkgs.mkShell {
    buildInputs = [
      rustPlatform.rust.cargo
      customPython
    ];

    shellHook = ''
       # workaround for https://github.com/NixOS/nixpkgs/blob/48dfc9fa97d762bce28cc8372a2dd3805d14c633/doc/languages-frameworks/python.section.md#python-setuppy-bdist_wheel-cannot-create-whl
       export SOURCE_DATE_EPOCH=315532800 # 1980
    '';
  }
