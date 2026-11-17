FoamFile { version 2.0; format ascii; class dictionary; object topoSetDict; }
actions
(
    { name Bcells; type cellSet; action new; source boxToCell; box (-4 -10 -15) (5 1 15); }
    { name SZone; type faceZoneSet; action new; source cellToFaceZone; set Bcells; }
);
