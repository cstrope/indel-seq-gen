LINEAGES = 
{
	Clade1, Clade2: "Ruminants" #mHKY,fnuc.freq# {5,0.02,idLD};
	Clade3: "Archaeoptyrix" #mF84# {0,0};
	Clade4: {7,0.0001/0.01, idLD/newLD} #p#;
	Taxon1(1,2): "Special_taxa"#i0.01# {2,0.0};
}

MOTIFS =
{
	Clade1:
		MARKER=a;
		NAME=Pribnow box;
		PATTERN=T-A-[TC]-{G}-[AT]-T;

	Clade3:
		MARKER=b;
		NAME= Removing constraint on TATA element.
		PATTERN=x(6)-x-x(0,3)-C-A-T;

	root:
		MARKER=c;
		NAME  = TATA box;
		PATTERN = T-A-T-A-A-[AT];
}
