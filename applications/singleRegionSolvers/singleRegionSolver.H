/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014 Johan Heyns - CSIR, South Africa
    Copyright (C) 2004-2011 OpenCFD Ltd.

-------------------------------------------------------------------------------
License
    This file is part of HiSA.

    HiSA is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HiSA is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with HiSA.  If not, see <http://www.gnu.org/licenses/>.

Application
    singleRegionSolver

Description
    Loads and calls the appropriate solver module for single region
    cases

See also
    multiRegionSolver

Authors
    Oliver Oxtoby
    Johan Heyns
    Council for Scientific and Industrial Research, South Africa

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pimpleControl.H"
#include "pseudotimeControl.H"
#include "solverModule.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    
    word regionName = fvMesh::defaultRegion;
    fvMesh mesh(
        IOobject(
            regionName,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ));

    dictionary dict;
    dict.add("solver", MODULE);
    autoPtr<solverModule> solverPtr =
        solverModule::New(Foam::polyMesh::defaultRegion, runTime, mesh, dict);

    solverModule& solver = solverPtr();

    solver.initialise();

    bool steadyState = solver.steadyState();

    if (!steadyState)
    {
        #include "createTimeControls.H"

        scalar scale = solver.timeStepScaling(maxCo);
        //Generate a pretend Courant number so we can use setInitialDeltaT.H unmodified
        scalar CoNum = maxCo/scale;
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    if (steadyState)
    {
        Info << "\nStarting pseudotime iteration loop\n" << endl;
    }
    else
    {
        Info<< "\nStarting time loop\n" << endl;
    }

    scalar prevCpuTime = runTime.elapsedCpuTime();
    while (runTime.run())
    {
        if (!steadyState)
        {
            #include "createTimeControls.H"

            //Find out what to scale deltaT by
            scalar scale = solver.timeStepScaling(maxCo);
            //Generate a pretend Courant number so we can use setDeltaT.H unmodified
            scalar CoNum = maxCo/scale;
            #include "setDeltaT.H"
        }

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        solver.beginTimeStep();

        // --- Outer corrector loop
        bool notFinished;
        #ifndef FOAM_VERSION_6
        while ( (notFinished = solver.solnControl().loop()) )
        {
        #else
        while (true)
        {
            // loop() function was removed from the base class
            if (isA<pimpleControl>(solver.solnControl()))
            {
                notFinished = refCast<pimpleControl>(solver.solnControl()).loop();
            }
            else
            {
                notFinished = refCast<pseudotimeControl>(solver.solnControl()).loop();
            }
            if (!notFinished)
            {
                break;
            }
        #endif
            solver.outerIteration();
            if (steadyState)
            {
                break;
            }
        }
        if (steadyState && !notFinished)
        {
            runTime.writeAndEnd();
        }
        else
        {
            runTime.write();
        }

        scalar cpuTime = runTime.elapsedCpuTime();
        Info<< "ExecutionTime = " << cpuTime << " s"
            << "  ExecutionStepTime = " << cpuTime-prevCpuTime << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        prevCpuTime = runTime.elapsedCpuTime();

    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
