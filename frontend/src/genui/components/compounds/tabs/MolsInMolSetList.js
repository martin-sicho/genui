import { Col, Row } from 'reactstrap';
import { CompoundListFromAPI } from '../../../index';
import React from 'react';

export default function MolsInMolSetList(props) {
  const molset = props.molset;
  const updateCondition = props.updateCondition ? props.updateCondition : (prevProps, nextProps) => {
    // console.log(prevProps.tasksRunning);
    // console.log(nextProps.tasksRunning);
    // console.log('xxx');
    return prevProps.tasksRunning !== nextProps.tasksRunning;
  };
  return (
    <Row>
      <Col sm="12">
        <h4>Molecules in {molset.name}</h4>
        <CompoundListFromAPI
          {...props}
          activitySetsIDs={molset.activities}
          showInfo={props.showInfo}
          updateCondition={updateCondition}
        />
      </Col>
    </Row>
  );
}