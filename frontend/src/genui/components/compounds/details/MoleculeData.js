import React from 'react';
import { Card, CardBody, CardHeader } from 'reactstrap';
import './compound-list-styles.css'

export function DataPair(props) {

  return (
    <p><b>{props.title}:</b> {props.data}</p>
  )
}

export function MoleculeData(props) {
  const mol = props.mol;
  console.log(mol);

  return (
    <React.Fragment>
      <Card className="compound-list-card">
        <CardHeader>
          <h3>Info</h3>
        </CardHeader>
        <CardBody>
          <DataPair title="SMILES" data={mol.smiles}/>
          <DataPair title="InChiKey" data={mol.inchiKey}/>
          {
            Object.keys(mol.extraArgs).map(key => (
              <DataPair key={key} title={key} data={mol.extraArgs[key].toString()} />
            ))
          }
        </CardBody>
      </Card>
    </React.Fragment>
  )
}