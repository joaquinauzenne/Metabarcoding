#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# scripts/lib/parse_yaml.sh
#
# Minimal, dependency-free YAML → shell-variable loader for the subset of YAML
# used by config/pipeline.yaml (nested scalar keys). Flattens nested keys with
# underscores under a prefix, e.g.:
#     blast:
#       perc_identity: 97      →   CFG_blast_perc_identity="97"
#
# Usage (source, then eval):
#     source "$(dirname "$0")/lib/parse_yaml.sh"
#     eval "$(parse_yaml path/to/config.yaml CFG_)"
#     echo "$CFG_blast_perc_identity"
#
# Scope / limitations (paired with capability, per project convention):
#   - DOES flatten nested scalar mappings and strip trailing "# ..." comments.
#   - Does NOT parse YAML sequences/lists, anchors, or multi-line scalars. The
#     `primers:` mapping is consumed by R (yaml::read_yaml), not by this loader.
#   - Keys must match [A-Za-z0-9_]; hyphenated keys (e.g. Riaz12S-V5) are skipped
#     here (harmless — bash never consumes them).
# ─────────────────────────────────────────────────────────────────────────────

parse_yaml() {
  local yaml_file="$1" prefix="${2:-}"
  local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs
  fs="$(printf '\034')"

  if [[ ! -f "$yaml_file" ]]; then
    printf "[ERROR] parse_yaml: config file not found: %s\n" "$yaml_file" >&2
    return 1
  fi

  sed -ne "s|^\(${s}\):|\1|" \
      -e "s|^\(${s}\)\(${w}\)${s}:${s}[\"']\(.*\)[\"']${s}\$|\1${fs}\2${fs}\3|p" \
      -e "s|^\(${s}\)\(${w}\)${s}:${s}\(.*\)${s}\$|\1${fs}\2${fs}\3|p" \
      "$yaml_file" |
  awk -F"$fs" -v prefix="$prefix" '{
    indent = length($1) / 2;
    vname[indent] = $2;
    for (i in vname) { if (i > indent) { delete vname[i] } }
    val = $3;
    # strip trailing inline "# comment" and surrounding whitespace
    sub(/[[:space:]]+#.*$/, "", val);
    sub(/[[:space:]]+$/, "", val);
    sub(/^[[:space:]]+/, "", val);
    if (length(val) > 0) {
      vn = "";
      for (i = 0; i < indent; i++) { vn = vn vname[i] "_" }
      gsub(/[^a-zA-Z0-9_]/, "_", vn);
      printf("%s%s%s=\"%s\"\n", prefix, vn, $2, val);
    }
  }'
}
