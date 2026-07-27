b10151203772: ~>> blockMesh
ls
cd 
ls
ls
blockMesh
checkMesh
snappyHexMesh -overwrite
exit
pwd
grep -n "mergeTolerance" system/snappyHexMeshDict
cp /Users/glebegor/Programming/projects/my/quantum-cfd/openfoam/system/snappyHexMeshDict system/snappyHexMeshDict
snappyHexMesh -overwrite
cd /home/openfoam
cat > system/snappyHexMeshDict <<'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      snappyHexMeshDict;
}

geometry
{
    wing_naca_0012.stl
    {
        type triSurfaceMesh;
        file "constant/triSurface/wing_naca_0012.stl";
    }
}

castellatedMeshControls
{
    maxLocalCells        100000;
    maxGlobalCells       2000000;
    mergeTolerance       1e-6;
    minRefinementCells   10;
    features
    (
    );

    refinementSurfaces
    {
        wing_naca_0012.stl
        {
            level           (2 2);
        }
    }

    resolveFeatureAngle  30;
    locationInMesh      (0.5 0 0.05);
    allowFreeStandingZoneFaces true;
}

snapControls
{
    nSmoothPatch         3;
    tolerance            2.0;
    nSolveIter           30;
    nRelaxIter           5;
}

addLayersControls
{
    relativeSizes       true;
    layers
    {
    }
    expansionRatio      1.0;
    finalLayerThickness 0.3;
    minThickness        0.1;
}

meshQualityControls
{
    maxNonOrtho 65;
    maxBoundarySkewness 20;
    maxInternalSkewness 4;
    maxConcave 80;
    minVol 1e-13;
    minTetQuality 1e-30;
    minArea -1;
    minTwist 0.02;
    minDeterminant 0.001;
    minFaceWeight 0.02;
    minVolRatio 0.01;
}

debug 0;
EOF

blockMesh
checkMesh
snappyHexMesh -overwrite
exit
ls
blockMesh
checkMesh
snappyHexMesh -overwrite
cd /home/openfoam
cat > system/snappyHexMeshDict <<'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      snappyHexMeshDict;
}

castellatedMesh true;
snap            true;
addLayers       false;

geometry
{
    wing_naca_0012.stl
    {
        type triSurfaceMesh;
        name wing;
    }
}

castellatedMeshControls
{
    maxLocalCells        100000;
    maxGlobalCells       2000000;
    minRefinementCells   10;
    nCellsBetweenLevels  2;

    features
    (
    );

    refinementSurfaces
    {
        wing
        {
            level (2 2);

            patchInfo
            {
                type wall;
            }
        }
    }

    resolveFeatureAngle 30;

    refinementRegions
    {
    }

    locationInMesh (5 0.5 0.055);

    allowFreeStandingZoneFaces true;
}

snapControls
{
    nSmoothPatch 3;
    tolerance    2.0;
    nSolveIter   30;
    nRelaxIter   5;

    nFeatureSnapIter    10;
    implicitFeatureSnap true;
    explicitFeatureSnap false;
}

addLayersControls
{
    relativeSizes true;

    layers
    {
    }

    expansionRatio      1.0;
    finalLayerThickness 0.3;
    minThickness        0.1;
}

meshQualityControls
{
    maxNonOrtho          65;
    maxBoundarySkewness  20;
    maxInternalSkewness  4;
    maxConcave           80;

    minVol               1e-13;
    minTetQuality        1e-30;
    minArea             -1;
    minTwist             0.02;
    minDeterminant       0.001;
    minFaceWeight        0.02;
    minVolRatio          0.01;

    nSmoothScale         4;
    errorReduction       0.75;
}

mergeTolerance 1e-6;

debug 0;
EOF

grep -n -B3 -A3 "mergeTolerance" system/snappyHexMeshDict
blockMesh
checkMesh
snappyHexMesh -overwrite
ls
snappyHexMesh -overwrite
blockMesh
checkMesh
snappyHexMesh -overwrite
blockMesh
checkMesh
snappyHexMesh -overwrite
blockMesh
checkMesh
snappyHexMesh -overwrite
blockMesh
checkMesh
snappyHexMesh -overwrite
pwd
ls
exit
