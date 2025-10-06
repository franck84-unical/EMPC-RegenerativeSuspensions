#!/usr/bin/env bash
set -euo pipefail

# 1) Create an empty GitHub repo named: EMPC-RegenerativeSuspensions
# 2) Replace USER below with your username or org
# 3) Run this script from the repo root (Codespaces or local)

git init
git checkout -b main
git lfs install || true

git add .
git commit -m "Initial commit: import EMPC regenerative suspensions package and paper"

git remote add origin git@github.com:USER/EMPC-RegenerativeSuspensions.git
# Or HTTPS:
# git remote add origin https://github.com/USER/EMPC-RegenerativeSuspensions.git

git push -u origin main
