import { Card, CardBody, CardImg } from 'reactstrap';
import React from 'react';
import './compound-list-styles.css'

export function MoleculePic(props) {
  const pic = props.mol.pics.find(pic => pic.image !== null);
  const As = props.as;

  const { as, ...rest } = props;
  return (
    pic ? <As {...rest} src={pic.image}/> : <p>No image found.</p>
  )
}

export function MoleculeDetail(props) {
  const mol = props.mol;

  return (
    <Card className="compound-list-card">
      <CardBody>
        <MoleculePic mol={mol} as={CardImg} top width="100%" alt={mol.smiles}/>
      </CardBody>
    </Card>
  )
}