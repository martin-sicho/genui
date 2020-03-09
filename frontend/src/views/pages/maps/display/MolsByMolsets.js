import React from 'react';
import {
  Col,
  Row,
} from 'reactstrap';
import {
  ActivitySetList,
  groupByMolset,
  MoleculeActivityDetail,
  MoleculeData,
  MoleculeDetail,
  TabWidget,
} from '../../../../genui';

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
            <Col md={3} sm={3}>
              <MoleculeData {...props} mol={mol}/>
            </Col>
            <Col md={6} sm={6}>
              <MoleculeActivityDetail
                {...props}
                mol={mol}
                component={ActivitySetList}
              />
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
          // points={data.map(item => item.point)}
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
  const molsets = props.molsets;
  const molsGroupedbyMolSet = groupByMolset(mols, molsets, (mol, idx) => {
    return {
      mol: mol
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