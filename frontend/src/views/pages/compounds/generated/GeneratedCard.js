import React from 'react';
import { GenericMolSetCard, MolSetTasks, LiveObject } from '../../../../genui';
import { Col, Row } from 'reactstrap';

function MolSetInfo(props) {
  const molset = props.molset;
  const molecules = props.molecules;

  return (
    <React.Fragment>
      <h4>Description</h4>
      <p>{molset.description}</p>

      <h4>Compounds</h4>
      <p>Unique in Total: {molecules.count}</p>
    </React.Fragment>
  )
}

function LiveMolSetInfo(props) {

  return (
    <LiveObject
      {...props}
      url={props.moleculesURL}
    >
      {
        molecules => {
          return <MolSetInfo {...props} molecules={molecules}/>
        }
      }
    </LiveObject>
  )
}

function GeneratedInfoTab(props) {

  return (
    <Row>
      <Col sm="12">
        <LiveMolSetInfo
          {...props}
        />
        <MolSetTasks
          progressURL={props.apiUrls.celeryProgress}
          tasks={props.tasks}
          molset={props.molset}
        />
      </Col>
    </Row>
  )
}

function GeneratedCard(props) {
  const tabs = [
    {
      title : "Info",
      renderedComponent : GeneratedInfoTab,
    }
  ];

  return (
    <GenericMolSetCard {...props} tabs={tabs}/>
  )
}

export default GeneratedCard;