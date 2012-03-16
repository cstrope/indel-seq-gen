LINEAGES = 
{
	Clade1, Clade2: "Ruminants" #mPAM,faaf2.freq# {5,0.02,idLD};
	Clade3: "Archaeoptyrix" #mPAM# {0,0};
	Clade4: {7,0.0001/0.01, idLD/newLD};
	Taxon1(1,2): "Special_taxa"#i0.01# {2,0.0};
	Taxon4: "Super_Special_Taxa"#mWAG#;
}

MOTIFS =
{
	Clade1:
		MARKER=a;
		NAME=Thioredoxin_1: Thioredoxin family active site, prosite pattern PS00194;
		PATTERN=[LIVMF]-[LIVMSTA]-x-[LIVMFYC]-[FYWSTHE]-x(2)-[FYWGTN]-C-[GATPLVE]-[PHYWSTA]-C-{I}-x-{A}-x(3)-[LIVMFYWT];

	Clade3:
		MARKER=b;
		NAME=G_PROTEIN_RECEP_F1_1,PS00237;G-protein coupled receptors family 1 signature;
		PATTERN=[GSTALIVMFYWC]-[GSTANCPDE]-{EDPKRH}-x-{PQ}-[LIVMNQGA]-{RK}-{RK}-[LIVMFT]-[GSTANC]-[LIVMFYWSTAC]-[DENH]-R-[FYWCSH]-{PE}-x-[LIVM];

	Taxon1:
		MARKER=c;
		NAME  = Thiol or something;
		PATTERN = C-x(2)-C;
}
