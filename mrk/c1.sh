## settings -------------------------------------------------------------------

### with mitochondrial sequence in one haplotype?
with_mito="FALSE"

### the gff files repository
dir_anno="${HOME}/data/nano-assemblies-pansn-2024/annotations"

## clmnt ----------------------------------------------------------------------

### folders
dir_base=$(dirname "${PWD}")
dir_out="${dir_base}/anno/gff"
if [[ -d "${dir_out}" ]]; then rm -rf "${dir_out}"; fi
mkdir -p "${dir_out}"

### the vector with the strain-haplotypes ids
declare -a gff_ids
file_ids="${dir_base}/ids/ids-ps.txt"
while IFS="\n" read -r one_line; do
    gff_ids+=("${one_line}")
done < "${file_ids}"

for strainhaplo_id in "${gff_ids[@]}"; do
  path_gff=$(find "${dir_anno}" -name "${strainhaplo_id}*gff")
  name_gff=$(basename "${path_gff}")
  strain_id=$(basename "${path_gff}" | cut -f 1 -d "-") # e.g. ADE
  mito_name="${strain_id}-mt-features.gff"
  path_mito_gff=$(find "${dir_anno}" -name "${mito_name}")
  gff_type=$(basename "${path_gff}" | cut -f 2 -d "-") # e.g. 0 or 1
  if [[ "${with_mito}" == "TRUE" && \
        -f "${path_mito_gff}" && \
        "${asse_type}" != "2" ]]; then
    cat <(grep -v "^#" "${path_gff}") <(grep -v "^#" "${path_mito_gff}") \
    > "${dir_out}/${name_gff}"
  else
    grep -v "^#" "${path_gff}" > "${dir_out}/${name_gff}"
  fi
done
