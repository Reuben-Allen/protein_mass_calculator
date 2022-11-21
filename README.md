# Protein Mass Calculator
This simple calculator will return the monoisoptopic and average mass of a peptide
along with the mass of the same protein enriched in a chosen isotope. An isotope can
be selected from a list of commonly utilized isotopic labels in biochemistry or a custom
isotope can be specied.  
To get started, simply run the code and follow the prompts.
## Example:
Let's find the mass of the protein SUMO-3 from *Homo sapiens* after enrichment with
nitrogen 15.
```
Enter the amino acid sequence as a continuous string of single letter abbreviations:
SEEKPKEGVKTENDHINLKVAGQDGSVVQFKIKRHTPLSKLMKAYCERQGLSMRQIRFRFDGQPINETDTPAQLEMEDEDTIDVFQQQTGG
Enter the atomic symbol of the desired isotope (multiple labels should be separate by ";"):
Key: "N" - Nitrogen 15
     "O" - Oxygen 18
     "H" - Deuterium
     "C" - Carbon 13
     "S" - Sulfur 34
     Custum: (H,C,O,N,S),MW
N

Molecular Formula: C-447 H-716 O-146 N-130 S-4
Average Mass (amu)
  Normal: 10393.5822
  Enriched: 10522.7255
Monoisotopic Mass (amu)
  Normal: 10387.1575
  Enriched: 10516.7727
```
