import React from 'react';
import {
  Card,
  CardBody,
  CardHeader,
  CardImg,
  Col,
  Row,
} from 'reactstrap';
import { groupByMolset, TabWidget } from '../../../../genui';
import './compound-list-styles.css'

function MoleculePic(props) {
  const pic = props.mol.pics.find(pic => pic.image !== null);
  const As = props.as;

  const { as, ...rest } = props;
  return (
    pic ? <As {...rest} src={pic.image}/> : <p>No image found.</p>
  )
}

function MoleculeDetail(props) {
  const mol = props.mol;

  return (
    <Card className="compound-list-card">
      <CardBody>
        <MoleculePic mol={mol} as={CardImg} top width="100%" alt={mol.smiles}/>
      </CardBody>
    </Card>
  )
}

function DataPair(props) {

  return (
    <p><b>{props.title}:</b> {props.data}</p>
  )
}

function MoleculeData(props) {
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

function MolSetDetail(props) {
  const mols = props.mols;

  return (
    <React.Fragment>
      {
        mols.map(mol => (
          <Row key={mol.id}>
            <Col md={3} sm={3}>
              <MoleculeDetail mol={mol}/>
            </Col>
            <Col md={9} sm={9}>
              <MoleculeData {...props} mol={mol}/>
            </Col>
          </Row>
        ))
      }
    </React.Fragment>
  )

}

export function MolsByMolsetsTabs(props) {
  const groupedMols = props.groupedMols;
  const molsets = props.molsets;
  const mols = props.mols;

  const tabs = [];
  let activeTab = undefined;
  Object.keys(groupedMols).forEach(molsetID => {
    const molset = molsets.find(item => item.id === Number(molsetID));
    const data = groupedMols[molsetID];
    if (!activeTab) {
      activeTab = molset.name;
    }
    tabs.push({
      title: molset.name
      , renderedComponent : props => (
        <MolSetDetail
          {...props}
          molset={molset}
          mols={data.map(item => item.mol)}
          points={data.map(item => item.point)}
        />
      )
    })
  });

  return (
    <div className="genui-map-molsets-grouped">
      {
        mols.length > 0 ? <TabWidget {...props} tabs={tabs} activeTab={activeTab}/> : <p>Select molecules in the map to see details.</p>
      }
    </div>
  )
}

export function MolsByMolsets(props) {
  const mols = props.selectedMols;
  const points = props.selectedPoints;
  const molsets = props.map.molsets;
  const molsGroupedbyMolSet = groupByMolset(mols, molsets, (mol, idx) => {
    return {
      mol: mol,
      point: points[idx]
    }
  });

  const RenderedComponent = props.component;
  return (
    <RenderedComponent
      {...props}
      molsets={molsets}
      mols={mols}
      groupedMols={molsGroupedbyMolSet}
    />
  )
}