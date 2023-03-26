sort $ids_to_download | uniq > ids_to_download_uniq.txt

for id in \$(cat ids_to_download_uniq.txt); do wget https://alphafold.ebi.ac.uk/files/\$id-model_v4.pdb; done

# Horrible coding - please change if keeping
# Generate links from name of the protein to the one matched_ref.pdb
python3 "${path_templates}/scripts/getlinks.py" ${template} "make_links_tmp.sh" "pdb"
[ -f ./make_links_tmp.sh ] && tr ', ' ' ' < make_links_tmp.sh > make_links.sh
[ -f ./make_links.sh ] && bash ./make_links.sh
[ -f ./make_links.sh ] && cat ./make_links.sh