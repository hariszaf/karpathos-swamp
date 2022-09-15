import sys

osmoadapt_kos_file  = open("kos-related-to-osmoadaptaion-only-term.tsv", "r")
mags_abund_taxonomy = open("mags_abundances_per_sample.csv", "r")
kos_per_mag_anvio   = open("kegg-metabolism-ko_hits-MATRIX.txt", "r")

number_of_osmo_kos_per_mag = {}
number_of_mags_per_ko      = {}
mag_ids                    = []
osmoadapt_kos              = []
mags_with_osmo_kos         = []

# Get MAG ids
counter = 0
for line in kos_per_mag_anvio:
    line = line.split("\t")
    if counter == 0:
        for mag_id in line[1:]:
            mag_ids.append(mag_id)
    else:
        break
    counter += 1

mag_ids[-1] = mag_ids[-1][:-1]


# Get osmoaddaptation-related KOs
for line in osmoadapt_kos_file: 
    osmoadapt_kos.append(line[:-1])

# Get MAGs with osmoadapt KOs
uniq_mag_ids_with_osmo_terms = set()
counter = 0
osmoadapt_kos_copy = osmoadapt_kos.copy()
for line in kos_per_mag_anvio:
    if counter > 1:
        line = line.split("\t")
        ko = line[0]
        if ko in osmoadapt_kos: 
            osmoadapt_kos_copy.remove(ko)
            line[-1] = line[-1][:-1]
            mag_id_counter = 0
            num_of_mags_with_the_ko = 0
            for abundance in line[1:]:
                if abundance != "0": 
                    num_of_mags_with_the_ko += 1
                    # print(ko, mag_ids[mag_id_counter] ,abundance)
                    uniq_mag_ids_with_osmo_terms.add(mag_ids[mag_id_counter])
                    if mag_ids[mag_id_counter] not in number_of_osmo_kos_per_mag:
                        number_of_osmo_kos_per_mag[mag_ids[mag_id_counter]] = 1
                    else:
                        number_of_osmo_kos_per_mag[mag_ids[mag_id_counter]] += 1                        
                mag_id_counter += 1
            number_of_mags_per_ko[ko] = num_of_mags_with_the_ko
    counter += 1

sorted_kos = sorted(number_of_osmo_kos_per_mag.items(), key=lambda x: x[1], reverse=True)

m_of_kos_in_n_mags = {}
for i in sorted_kos:
    if i[1] not in m_of_kos_in_n_mags:
        m_of_kos_in_n_mags[i[1]] = 1
    else:
        m_of_kos_in_n_mags[i[1]] += 1

for i, j in m_of_kos_in_n_mags.items():
    print(i, "\t", j)

# print(osmoadapt_kos_copy)
