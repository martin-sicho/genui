import React from 'react';
import './compound-list-styles.css'

export function DataPair(props) {

  return (
    <p><b>{props.title}:</b> {props.data}</p>
  )
}

export function MoleculeMetadata(props) {
  const mol = props.mol;

  return (
    <React.Fragment>
      <DataPair title="SMILES" data={mol.smiles}/>
      <DataPair title="InChiKey" data={mol.inchiKey}/>
      {
        Object.keys(mol.extraArgs).map(key => (
          <DataPair key={key} title={key} data={mol.extraArgs[key].toString()} />
        ))
      }
    </React.Fragment>
  )
}