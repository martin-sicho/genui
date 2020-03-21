import React from 'react';
import { Col, Row } from 'reactstrap';
import { ActivitiesAggregator, groupBy, MolSetsTabs, MolsToMolSetGroups, TabWidget } from '../../../../genui';
import ActivitySummaryPlotter from './ActivitySummaryPlotter';

export default function SelectedPage(props) {
  const selectedMols = props.selectedMols;
  return (
    selectedMols.length > 0 ? (<React.Fragment>
        <h1>Selected Compounds in Map</h1>
        <hr/>

        <Row>
          <Col sm={12}>
            <h2>Activity Summary</h2>
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
                        <ActivitySummaryPlotter
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
          </Col>
        </Row>

        <hr/>

        <Row>
          <Col sm={12}>
            <h2>List of Selected Compounds</h2>
            <MolsToMolSetGroups
              {...props}
              mols={selectedMols}
            >
              {
                (groups) => {
                  return (
                    <MolSetsTabs
                      {...props}
                      groupedMols={groups}
                      mols={selectedMols}
                    />
                  )
                }
              }
            </MolsToMolSetGroups>
          </Col>
        </Row>
      </React.Fragment>
    ) : <p>No compounds are selected, yet. You can select compounds in the map with the lasso tool, for example.</p>
  )
}