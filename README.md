# The Alice Assembler

The Alice assembler is named after the character of Lewis Caroll's _Alice in Wonderland_, and more precisely to the "drink-me potion" and "eat-me cake", which make Alice respectively shrink and grow. The idea of the Alice assembler is to shrink the input reads, perform plain De Bruijn Graph assembly on shrunken data and inflate the obtained assembly back to normal size. The compression is a Mapping-friendly Sequence Reduction (MSR) of high order, which allows all base to be taken into account during compression.
![alice_compression](https://github.com/rolandfaure/Alice-asm/blob/master/alice_compression.png)
