{
  description = "DMS-nix";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/b134951a4c9f3c995fd7be05f3243f8ecd65d798"; # nixpkgs = nixos-24.05
    rpkgs.url = "github:rstats-on-nix/nixpkgs/c2c46d90ff2d7e7f2c4982cab38de3bb3cab0225"; # rstats-on-nix = 2025-01-14
    cmdstanr.url = "github:stan-dev/cmdstanr/02259ef7aa2a8b1c8de2fa3fc42a9feafd789288"; # cmdstanr = v0.8.1
    cmdstanr.flake = false;
  };

  outputs = { self, nixpkgs, rpkgs, cmdstanr }:
    
    let
      
      system = "x86_64-linux";
      pkgs = nixpkgs.legacyPackages.${system};
      rpkgs_pkgs = rpkgs.legacyPackages.${system};
      
      r_packages = builtins.attrValues {
        inherit (rpkgs_pkgs.rPackages) 
          argparse
          bench
          brms
          broom
          broom_mixed
          checkmate
          data_table
          drc
          emmeans
          factoextra
          furrr
          future
          future_callr
          ggbeeswarm
          ggh4x
          ggrastr
          ggrepel
          glmmTMB
          IRkernel
          janitor
          jsonlite
          Matrix
          patchwork
          posterior
          processx
          R6
          rix
          scales
          scico
          tidyverse
          TMB
          withr;
      };

      system_packages = builtins.attrValues {
        inherit (rpkgs_pkgs) 
          nix
          R;
        inherit (pkgs)
          cmdstan;
      };

      r_cmdstanr = pkgs.rPackages.buildRPackage {
        name = "cmdstanr";
        src = cmdstanr;
        propagatedBuildInputs = r_packages;
      };

      python_packages = builtins.attrValues {
        inherit (pkgs.python312Packages)
          jupyter-core
          jupyterlab;
      };

    in

      {
        devShells.x86_64-linux.default = pkgs.mkShell {
          buildInputs = [ r_packages system_packages python_packages r_cmdstanr ];
          shellHook = ''
            export CMDSTAN="${pkgs.cmdstan}/opt/cmdstan/"
          '';
        };
      };
}

