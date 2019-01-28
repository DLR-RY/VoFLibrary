/*---------------------------------------------------------------------------*\
            Copyright (c) 2017-2019, German Aerospace Center (DLR)
-------------------------------------------------------------------------------
License
    This file is part of the VoFLibrary source code library, which is an 
	unofficial extension to OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::zoneCPCStencil

Description
    computes a cell point cell stencil in a narrow band. resizes in case
    of topological change

Author
    Henning Scheufler, DLR, all rights reserved.
    
SourceFiles
    zoneCPCStencil.C

\*---------------------------------------------------------------------------*/

#ifndef zoneCPCStencil_H
#define zoneCPCStencil_H

#include "zoneCellStencils.H"
#include "boolList.H"
#include "HashSet.H"
#include "Map.H"
#include "indirectPrimitivePatch.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class zoneCPCStencil Declaration
\*---------------------------------------------------------------------------*/

class zoneCPCStencil
:
    public zoneCellStencils
{
    // Private Member Functions
        labelList nonEmptyBoundaryPoints_;

        boolList isValidBFace_;

        boolList uptodate_;

        Map<bool> syncCouledBoundaryPoints
        (
            const boolList& zone,
            const labelList& boundaryPoints
        ) const;

        //- Calculates per point the neighbour data (= pointCells)
        void calcPointBoundaryData
        (
            const boolList& zone,
            const boolList& isValidBFace,
            const labelList& boundaryPoints,
            Map<labelList>& neiGlobal
        ) const;

        void calculateStencil
        (
            const boolList& zone,
            labelListList& globalCellCells
        );



public:

    // Constructors
        // Declare name of the class and its debug switch
        TypeName("zoneCPCStencil");

        //- Construct from all cells and boundary faces
        explicit zoneCPCStencil(const fvMesh&);

        //- Calculates per cell the neighbour data
        //  (= cell or boundary in global numbering).
        //  First element is always cell itself!

        virtual void updateMesh(const mapPolyMesh& mpm);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
