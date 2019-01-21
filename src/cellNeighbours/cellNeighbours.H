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
    Foam::cellNeighbours

Description
    Class that provides cell-point-cell stencil. A labelList<DynamicList <label>>
    is returned. The elements are in globalCellNumbering (cell and boundaryFaces)
    (see. class globalIndex). The boundary values and processor information are
    provided by a Map (used in plicRDF isoRDFIter)

Author
    Henning Scheufler, DLR, all rights reserved.

SourceFiles
    cellNeighbours.C
    cellNeighboursI.H

\*---------------------------------------------------------------------------*/

#ifndef cellNeighbours_H
#define cellNeighbours_H

#include "OFstream.H"
#include "cutCellIso.H"
#include "fvMesh.H"
#include "globalIndex.H"
#include "scalarMatrices.H"
#include "surfaceFields.H"
#include "token.H"
#include "volFields.H"

#include "leastSquareCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cellNeighbours Declaration
\*---------------------------------------------------------------------------*/

class cellNeighbours
{

  private:
    //- Combine operator for labelLists
    template <class type> class appendOp
    {
      public:
        void operator()(List<type>& x, const List<type>& y) const
        {
            if (y.size())
            {
                if (x.empty())
                {
                    x = y;
                }
                else
                {
                    x.append(y);
                }
            }
        }
    };

    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        // cell-point-cell stencil elements are in global addressing
        List<DynamicList<label>> stencil_;

        // parallel

        // point that are on a processorpatch
        boolList onProcPatch_;

        // Map and List providing information on which procssor a point is
        Map<label> procPointsMap_;
        DynamicList<labelPairList> procPoints_;

        // Map for getting unique indices
        labelHashSet cellSet_;

        // Map is synced to find points and other processor
        Map<bool> syncInterfacePointMap_;

        // Map provides global Cells Numbering (on points)
        Map<List<label>> syncNeiCellMap_;

        // Provides information of processor addressing
        Map<List<label>> boundaryCellAddressing_;

        // provides global addressing
        const globalIndex globalNumbering_;

        // appends new processor number to procs
        void appendUniqueProcs(const label& pI, DynamicList<label>& procs);

        // calcs cell stencil and updates boundaryCellAddressing_
        void calcLocalCellStencil(const label& cellI);

        // check if cellI is next to a processor and calcs boundaryCellAddressing_
        void calcBoundaryProcAddressing(const label& cellI);

         //- gives patchNumber and patchFaceNumber for a given Geometric volume field
        void patchIFace
        (
            const label faceI,
            label& patchNumber,
            label& patchFaceNumber
        ) const;


        //- gives patchNumber and patchFaceNumber for a given Geometric volume field
        template <typename Type>
        void faceValue
        (
            const GeometricField<Type, fvPatchField, volMesh>& f,
            const label facei,
            label& patchNumber,
            label& patchFaceNumber
        ) const;

        // syncs syncInterfacePointMap_ on processors
        void syncInterfacePoints(const boolList& nextToInterface);

        // gets adjacent global cell and updates syncNeiCellMap_
        void syncGlobalIdxMap(const boolList& nextToInterface);

        // marks on processor patches updates onProcPatch_
        void markPointPatches();

  public:
    // Constructors

        //- Construct from fvMesh
        cellNeighbours(const fvMesh&);

        // Member functions

        //- updates the stencil with boolList the size has to match
        // mesh.nCells()
        void updateStencil(const boolList& nextToInterface);

        //- Returns stenicl and provides a Map with globalNumbering
        template <typename Type>
        const List<DynamicList<label>>& getDataFromOtherProcs
        (
            const GeometricField<Type, fvPatchField, volMesh>& f,
            const boolList& nextToInterface,
            Map<Type>& map
        );

        template <typename Type>
        void getPointDataFromOtherProcs
        (
            const GeometricField<Type, fvPatchField, volMesh>& f,
            const boolList& interfacePoints,
            Map<List<Type>>& map
        );

        //- Returns stenicl and provides a Map with globalNumbering for a List of
        //geometricFields
        template <typename Type>
        const List<DynamicList<label>>& getDataFromOtherProcs
        (
            const PtrList<GeometricField<Type, fvPatchField, volMesh>>& f,
            const boolList& nextToInterface,
            Map<List<Type>>& map
        );


        //- return stencil reference
        const List<DynamicList<label>>& getStencil()
        {
            return stencil_;
        }

        //- return addressing reference
        const globalIndex& globalNumbering()
        {
            return globalNumbering_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "cellNeighboursI.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
