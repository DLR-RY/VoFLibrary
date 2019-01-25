/*---------------------------------------------------------------------------*\
    Modified work | Copyright (c) 2017-2019, German Aerospace Center (DLR)
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
    Foam::cutCellIso

Description
    Class for cutting a cell, celli, of an fvMesh, mesh_, at its intersection
    with an isosurface defined by the mesh point values f_ and the cutValue,
    isoValue_.

    Reference:
    \verbatim
        Roenby, J., Bredmose, H. and Jasak, H. (2016).
        A computational method for sharp interface advection
        Royal Society Open Science, 3
        doi 10.1098/rsos.160405

        Henning Scheufler, Johan Roenby,
        Accurate and efficient surface reconstruction from volume
        fraction data on general meshes,
        Journal of Computational Physics, 2019,
        doi 10.1016/j.jcp.2019.01.009

    \endverbatim
    
Author
    Johan Roenby, DHI, all rights reserved.
    Modifications: Henning Scheufler, DLR, all rights reserved.

SourceFiles
    cutCellIso.C

\*---------------------------------------------------------------------------*/

#ifndef cutCellIso_H
#define cutCellIso_H

#include "OFstream.H"
#include "cutCell.H"
#include "cutFaceIso.H"
#include "fvMesh.H"
#include "scalarMatrices.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cutCellIso Declaration
\*---------------------------------------------------------------------------*/

class cutCellIso
:
    public cutCell
{
    // Private data

        //- Mesh whose cells and faces to cut at their intersection with an
        //  isosurface.
        const fvMesh& mesh_;

        //- Cell to cut
        label cellI_;

        //- Isofunction values at mesh points. f_size() = mesh_.nPoints().
        scalarField& f_;

        //- Isovalue used to cut cell
        scalar cutValue_;

        //- An cutFaceIso object to get access to its face cutting functionality
        cutFaceIso cutFace_;

        //- List of face centres for CutFaces
        DynamicList<point> cutFaceCentres_;

        //- List of face area vectors for isoCutFaces
        DynamicList<vector> cutFaceAreas_;

        //- Storage for subFace edges belonging to isoFace
        DynamicList<DynamicList<point>> isoFaceEdges_;

        //- Points constituting the cell-isosurface intersection (isoface)
        DynamicList<point> facePoints_;

        //- Face centre of the cutFace
        point faceCentre_;

        //- Face normal of the isoface by convention pointing from high to low
        //  values (i.e. opposite of the gradient vector).
        vector faceArea_;

        //- Cell centre of the subcell of celli which is "fully submerged", i.e.
        //  where the function value is higher than the isoValue_
        point subCellCentre_;

        //- Volume of fully submerged subcell
        scalar subCellVolume_;

        //- Volume of Fluid for cellI (subCellVolume_/mesh_.V()[cellI])
        scalar VOF_;

        //- A cell status label taking one of the values:
        //
        //  - -1: cell is fully below the isosurface
        //  -  0: cell is cut
        //  - +1: cell is fully above the isosurface
        label cellStatus_;

    // Private Member Functions

  public:
    // Constructors

        //- Construct from fvMesh and a scalarField
        //  Length of scalarField should equal number of mesh points
        cutCellIso(const fvMesh&, scalarField& f);

        // Member functions

        //- sets internal values and returns face status
        label calcSubCell(const label& cellI, const scalar cutValue);

        //- returns subCellCentre
        point subCellCentre();

        //- returns subCellVolume
        scalar subCellVolume();

        //- returns the points of the cutting isoface
        DynamicList<point> facePoints();

        //- returns the centre of the cutting isoface
        point faceCentre();

        //- returns the area normal vector of the cutting isoface
        vector faceArea();

        //- returns cellStatus
        label cellStatus();

        //- returns volume of fluid value
        scalar VolumeOfFluid();

        //- returns cutValue
        scalar cutValue() const;

        //- resets internal values
        void clearStorage();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
