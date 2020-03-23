import { Card, CardBody, CardHeader, Col, Row } from 'reactstrap';
import { ActivitiesTable, MoleculeActivityProvider, MoleculeImage } from '../../..';
import React from 'react';

function ActivitySetFlatView(props) {
  const activities = props.activities;
  const actsets = props.activitySets;

  return (
    <React.Fragment>
      {
        Object.keys(actsets).map(actsetKey => {
          const actset = actsets[actsetKey];
          const actsetActivities = activities[actset.id];
          if (actsetActivities &&  actsetActivities.length > 0) {
            // const byType = groupBy(actsetActivities, 'type.id');
            // console.log(byType);
            actsetActivities.sort((item) => item.type.name);

            // FIXME: filtering should not be necessary, fix ComponentWithPagedResources or MoleculeActivityProvider so that it does not leak previous information
            const filteredActivities = [];
            actsetActivities.forEach(activity => {
              if (activity.molecule === props.mol.id) {
                filteredActivities.push(activity);
              }
            });
            const extraData = [
              {
                header: "Type",
                data: filteredActivities.map(activity => activity.type.value)
              }
            ];
            return (
              <Card key={actsetKey}>
                <CardHeader>{actset.name}</CardHeader>
                <CardBody>
                  <ActivitiesTable
                    activities={filteredActivities}
                    extraData={extraData}
                    extraDataAppend={false}
                  />
                </CardBody>
              </Card>
            )
          } else {
            return null;
          }
        })
      }
    </React.Fragment>
  )
}

export default function CompoundOverview(props) {
  const mol = props.mol;

  return (
    <React.Fragment>
      <Row>
        <Col sm={12}>
          <MoleculeImage
            {...props}
          />
        </Col>
      </Row>

      <hr/>

      <h3>Activities</h3>
      <Row>
        <Col sm={12}>
          <MoleculeActivityProvider
            {...props}
            mol={mol}
            updateCondition={(prevProps, currentProps) => {
              return prevProps.mol && (prevProps.mol.id !== currentProps.mol.id)
            }}
            component={ActivitySetFlatView}
          />
        </Col>
      </Row>

      <hr/>

      <h3>PhysChem Properties</h3>
      <Row>
        <Col sm={12}>
          <div>Fetch them here...</div>
        </Col>
      </Row>
    </React.Fragment>
  )
}