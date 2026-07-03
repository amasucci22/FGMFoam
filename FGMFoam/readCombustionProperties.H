Info<< "Reading combustion properties\n" << endl;

IOdictionary combustionProperties
(
	IOobject
	(
		"combustionProperties",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ_IF_MODIFIED,
		IOobject::NO_WRITE
	)
);


const word FGMModel //Run-time not allowed
(
	combustionProperties.lookup("FGMModel")
);

word PDModel //Run-time modifiable
(
	combustionProperties.lookup("PDModel")
);

Switch reacting //Run time modifiable 
(
	combustionProperties.lookup("ignite")
);

Switch betacDynamic //Run-time modifiable
(
	combustionProperties.lookup("betacDynamic")
);

dimensionedScalar betacVal //Dynamic computation will use a different variable; run-time modifiable
(
	combustionProperties.lookup("betaC")
);




    
// **************************************************************** //


