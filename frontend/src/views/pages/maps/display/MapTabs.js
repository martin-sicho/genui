import React from 'react';
import MapPage from './MapPage';
import { Col, Row } from 'reactstrap';
import { MolsToMolSetGroups, MolSetsTabs, TabWidget, groupBy, ActivitiesAggregator } from '../../../../genui';
import ActivitySummaryPlotter from './ActivitySummaryPlotter';

function SelectedPage(props) {
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

export default function MapTabs(props) {
  const [selectedMolsInMap, setSelectedMolsInMap] = React.useState([]);
  const [selectedMolsRevision, setSelectedMolsInMapRevision] = React.useState(0);
  const [selectedPointsInMap, setSelectedPointsInMap] = React.useState([]);
  // const [activeTab, setActiveTab] = React.useState('Map');

  // const tabs = {
  //   map: {
  //     title: "Map"
  //   },
  //   selected: {
  //     title: "Selected"
  //   }
  // };

  const tabs = [
    {
      title: "Map",
      renderedComponent: MapPage
    },
    {
      title: "Selection",
      renderedComponent: SelectedPage
    }
  ];

  return (
    <TabWidget
      {...props}
      tabs={tabs}
      activeTab={tabs[0].title}
      selectedMols={selectedMolsInMap}
      selectedPoints={selectedPointsInMap}
      onMolsSelect={setSelectedMolsInMap}
      selectedMolsRevision={selectedMolsRevision}
      setSelectedMolsInMapRevision={setSelectedMolsInMapRevision}
      onPointsSelect={setSelectedPointsInMap}
    />
  );

  // return (
  //   <Tabs
  //     defaultActiveKey={tabs.map.title}
  //     activeKey={activeTab}
  //     onSelect={k => setActiveTab(k)}
  //     id="noanim-tab-example"
  //   >
  //     <Tab
  //       eventKey={tabs.map.title}
  //       title={tabs.map.title}
  //       // className={classnames({ active: activeTab === summaries.map.title })}
  //     >
  //       <MapPage
  //         {...props}
  //         selectedMols={selectedMolsInMap}
  //         selectedPoints={selectedPointsInMap}
  //         onMolsSelect={setSelectedMolsInMap}
  //         selectedMolsRevision={selectedMolsRevision}
  //         setSelectedMolsInMapRevision={setSelectedMolsInMapRevision}
  //         onPointsSelect={setSelectedPointsInMap}
  //       />
  //     </Tab>
  //
  //     <Tab
  //       eventKey={tabs.selected.title}
  //       title={tabs.selected.title}
  //       // className={classnames({ active: activeTab === summaries.selected.title })}
  //     >
  //       <Selected
  //         {...props}
  //         selectedMols={selectedMolsInMap}
  //         selectedPoints={selectedPointsInMap}
  //         onMolsSelect={setSelectedMolsInMap}
  //         selectedMolsRevision={selectedMolsRevision}
  //         setSelectedMolsInMapRevision={setSelectedMolsInMapRevision}
  //         onPointsSelect={setSelectedPointsInMap}
  //       />
  //     </Tab>
  //
  //   </Tabs>
  // )
}