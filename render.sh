#!/usr/bin/env bash
# Re-run the full study and regenerate every output.
# HTML is rendered first because its self-contained mode deletes the shared
# figure folder; the GitHub-markdown render then recreates the figures it needs.
set -euo pipefail
cd "$(dirname "$0")/report"
quarto render clt_project_1.qmd --to html
quarto render clt_project_1.qmd --to gfm
