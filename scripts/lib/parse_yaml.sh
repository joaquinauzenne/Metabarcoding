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
#   - Keys must match [A-Za-z0-9_]. A key containing anything else (e.g. the
#     primer set Riaz12S-V5) cannot be represented as a shell variable, so it is
#     SKIPPED — along with its children. Without that, the unmatched parent line
#     was dropped while its children were still emitted, silently reattributing
#     them to the previous sibling key (Riaz12S-V5's min_len would overwrite
#     Taberlat's). Such keys are only consumed by R, which parses YAML properly.
# ─────────────────────────────────────────────────────────────────────────────

parse_yaml() {
  local yaml_file="$1" prefix="${2:-}"
  # `w` deliberately accepts MORE than a shell identifier (hyphens, dots) so that
  # keys like `Riaz12S-V5:` are still matched here and reach awk, which then skips
  # them and their children. If they were never matched, the parent line would be
  # dropped while its children still emitted under the previous sibling's name.
  local s='[[:space:]]*' w='[A-Za-z0-9_.-]*' fs
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
    key = $2;
    val = $3;

    # Leaving a nesting level invalidates every deeper level of BOTH maps.
    for (i in vname) { if ((i+0) > indent) { delete vname[i] } }
    for (i in bad)   { if ((i+0) > indent) { delete bad[i]   } }

    # A key that is not a valid shell identifier cannot be emitted, and neither
    # can anything nested under it -- mark this level bad and drop the line.
    if (key !~ /^[A-Za-z0-9_]+$/) { bad[indent] = 1; delete vname[indent]; next }
    delete bad[indent];
    vname[indent] = key;

    # Suppress this line if any ANCESTOR level was skipped, so orphaned children
    # are never reattributed to the previous sibling key.
    for (i = 0; i < indent; i++) { if (i in bad) next }

    # strip trailing inline "# comment" and surrounding whitespace
    sub(/[[:space:]]+#.*$/, "", val);
    sub(/[[:space:]]+$/, "", val);
    sub(/^[[:space:]]+/, "", val);
    if (length(val) > 0) {
      vn = "";
      for (i = 0; i < indent; i++) { vn = vn vname[i] "_" }
      gsub(/[^a-zA-Z0-9_]/, "_", vn);
      printf("%s%s%s=\"%s\"\n", prefix, vn, key, val);
    }
  }'
}
