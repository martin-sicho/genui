import { Col, Row } from 'reactstrap';
import { LiveObject, MolSetTasks } from '../../../../../genui';
import React from 'react';

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

function GeneratedInfo(props) {

  return (
    <Row>
      <Col sm="12">
        <LiveMolSetInfo
          {...props}
        />
        <MolSetTasks
          {...props}
          progressURL={props.apiUrls.celeryProgress}
        />
      </Col>
    </Row>
  )
}

export default GeneratedInfo;