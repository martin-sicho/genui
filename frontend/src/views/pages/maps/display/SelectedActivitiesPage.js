import React from 'react';
import { Col, Row } from 'reactstrap';
import { ActivitiesAggregator, groupBy, TabWidget } from '../../../../genui';
import ActivitySummary from './ActivitySummary';

export default function SelectedActivitiesPage(props) {
  const selectedMols = props.selectedMols;
  return (
     <React.Fragment>
        <h1>Activity Summary for Selected Compounds in the Map</h1>
        <hr/>

        <Row>
          <Col sm={12}>
            {
              selectedMols.length > 0 ? (
                <ActivitiesAggregator
                  {...props}
                  mols={selectedMols}
                  resourceUpdateCondition={(prevProps, currentProps) => prevProps.selectedMolsRevision !== currentProps.selectedMolsRevision}
                >
                  {
                    (activities) => {
                      if (activities) {
                        const groupedActivities = groupBy(activities, 'type.id');
                        // console.log(groupedActivities);

                        const tabs = groupedActivities.map(group => ({
                          title: group[0].type.value,
                          renderedComponent: (props) => (
                            <ActivitySummary
                              {...props}
                              type={group[0].type}
                              activities={group}
                            />
                          )
                        }));

                        return (
                          <TabWidget
                            {...props}
                            tabs={tabs}
                          />
                        )
                      } else {
                        return <div>Loading...</div>;
                      }
                    }
                  }
                </ActivitiesAggregator>
              ) : <p>Select compounds in the map to see details.</p>
            }
          </Col>
        </Row>
      </React.Fragment>
  )
}