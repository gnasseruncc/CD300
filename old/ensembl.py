import requests, sys

def get_species_data(species_name):
    server = "https://rest.ensembl.org"
    ext = f"/info/genomes/{species_name}?"
    
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    
    # if not r.ok:
    #   r.raise_for_status()
    
    decoded = r.json()
    return decoded

def main():
    species_list = [
        "Cynocephalus volans",
        "Galeopterus variegatus",
        "Tupaia glis",
        "Tupaia longipes",
        "Tupaia tana",
        "Dendrogale murina frenata",
        "Ptilocercus lowii",
        "Tupaia gracilis gracilis",
        "Tupaia minor",
        "Abrothrix_longipilis",
        "Ammospermophilus leucurus",
        "Aplodontia rufa",
        "Apodemus sylvaticus",
        "Cynomys ludovicianus",
        "Cavia porcellus",
        "Galea musteloides",
        "Tamias sibiricus",
        "Lariscus insignis",
        "Liomys pictus",
        "Marmota caligata",
        "Marmota himalayana",
        "Marmota monax",
        "Menetes berdmorei",
        "Mus musculus",
        "Neotoma bryanti",
        "Nesomys rufus",
        "Pedetes capensis",
        "Perognathus amplus",
        "Rattus rattus",
        "Sciurotamias davidianus",
        "Urocitellus richardsonii",
        "Tamias amoenus",
        "Xerus erythropus",
        "Xerus inauris",
        "Callosciurus erythraeus",
        "Chiropodomys gliroides",
        "Coendou prehensilis",
        "Dactylomys dactylinus",
        "Dendromus mystacalis",
        "Dremomys pernyi",
        "Eliurus webbi",
        "Erethizon dorsatum",
        "Funisciurus pyrropus",
        "Heliosciurus rufobrachium",
        "Lophiomys imhausi",
        "Micromys minutus",
        "Microsciurus alfari",
        "Microsciurus flaviventer",
        "Myosciurus pumilio",
        "Nannosciurus melanotis",
        "Paraxerus cepapi",
        "Paraxerus ochraceus",
        "Phloeomys cumingi",
        "Ratufa bicolor",
        "Sciurus aberti",
        "Sciurus carolinensis",
        "Sciurus niger",
        "Sundasciurus lowii",
        "Tamiasciurus douglasii",
        "Tamiasciurus hudsonicus",
        "Tamiops swinhoei",
        "Aeromys tephromelas",
        "Anomalurus beecrofti",
        "Eoglaucomys fimbriatus",
        "Glaucomys sabrinus",
        "Petaurista elegans",
        "Petaurista petaurista",
        "Petaurista philippensis",
        "Trogopterus xanthipes",
        "Aotus trivirgatus",
        "Arctocebus calabarensis",
        "Avahi laniger",
        "Plecturocebus moloch",
        "Callicebus personatus",
        "Callimico goeldii",
        "Callithrix geoffroyi",
        "Callithrix pygmaea",
        "Callithrix jacchus",
        "Cebus capucinus",
        "Cheirogaleus major",
        "Cheirogaleus medius",
        "Daubentonia madagascariensis",
        "Eulemur coronatus",
        "Eulemur fulvus",
        "Eulemur macaco",
        "Eulemur mongoz",
        "Eulemur rubriventer",
        "Euoticus elegantulus",
        "Galago alleni",
        "Galago moholi",
        "Galago senegalensis",
        "Galagoides zanzibaricus",
        "Galagoides demidoff",
        "Hapalemur griseus",
        "Prolemur simus",
        "Indri indri",
        "Lemur catta",
        "Leontopithecus chrysomelas",
        "Leontopithecus rosalia",
        "Lepilemur leucopus",
        "Lepilemur mustelinus",
        "Lepilemur ruficaudatus",
        "Loris lydekkerianus nordicus",
        "Loris tardigradus",
        "Microcebus murinus",
        "Mirza coquereli",
        "Nycticebus coucang",
        "Nycticebus pygmaeus",
        "Otolemur crassicaudatus",
        "Otolemur garnettii",
        "Perodicticus potto",
        "Phaner pallescens",
        "Propithecus diadema",
        "Propithecus verreauxi",
        "Saguinus midas",
        "Saguinus oedipus",
        "Saimiri sciureus",
        "Tarsius bancanus",
        "Tarsius spectrum",
        "Carlito syrichta",
        "Varecia variegata"
    ]

    for species in species_list:
        # Format species name for URL
        species_name1 = species.lower()
        species_name2 = species_name1.replace(" ", "_")
        data = get_species_data(species_name2)
        
        if data:
            display_name = data.get('display_name', 'N/A')
            scientific_name = data.get('scientific_name', 'N/A')
            assembly_accession = data.get('assembly_accession', 'N/A')
            print(f"Species: {species}, Display Name: {display_name}, Scientific Name: {scientific_name}, Assembly Accession: {assembly_accession}")

if __name__ == "__main__":
    main()
