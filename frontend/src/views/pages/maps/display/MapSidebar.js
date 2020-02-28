import React from 'react';

function MoleculeDetail(props) {
  const mol = props.mol;
  const point = props.point;
  const molsetColor = point.fullData.marker.color;

  return (
    <p>{mol.smiles}</p>
  )
}

export default function MapSidebar(props) {
  const mols = props.selectedMols;
  const points = props.selectedPoints;

  return (
    <div className="genui-map-sidebar">
      <h3>Selected Molecules</h3>

      {
        mols.length > 0 ? mols.map((mol, index) => {
          return <MoleculeDetail key={mol.id} mol={mol} point={points[index]}/>
        }) : <p>No molecules selected. Select them in the map to see details.</p>
      }
    </div>
  )
}