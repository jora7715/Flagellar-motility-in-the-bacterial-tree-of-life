####SIMULATE METAGENOMES FROM GENOME MIXES
rm(list=ls()) 

library(Biostrings)
library(ballgown)

simulate_metagenomes <- function(genome_file, fragment_length, num_fragments) {
  # Read the genome sequences
  genomes <- readDNAStringSet(genome_file)
  
  simulated_fragments <- DNAStringSet()
  for (i in 1:num_fragments) {
    # Choose a random genome
    genome_id <- sample(names(genomes))[sample(1:length(names(genomes)), 1)]
    genome <- genomes[genome_id]
    
    # Randomly choose a fragment start position
    start <- sample(1:(length(genome) - fragment_length), 1)
    
    # Extract the fragment
    fragment <- genome[start:(start + fragment_length - 1)]
    names(fragment) <- paste0(genome_id, "_", start, "_", start + fragment_length - 1)
    
    simulated_fragments <- c(simulated_fragments, fragment)
  }
  
  return(simulated_fragments)
}

# Example usage
genome_file <- "/data/ramonedaj/Chemotaxis/SimulatedMetagenomes/MetagenomeAnnotationExample/Allmotile/combined_allmotile.fna"
fragment_length <- 150
num_fragments <- 1000000

fragments <- simulate_metagenomes(genome_file, fragment_length, num_fragments)
writeXStringSet(fragments, file = "/data/ramonedaj/Chemotaxis/SimulatedMetagenomes/MetagenomeAnnotationExample/Allmotile/simulated_allmotile_metagenome.fasta")
