import React from 'react';
import { ComponentWithObjects, ComponentWithResources } from '../../../../genui';
import MapPage from './MapPage';
import { Col, Row } from 'reactstrap';
import { MolsByMolsets, MolsByMolsetsTabs } from './MolsByMolsets';
import { Tab, Tabs } from 'react-bootstrap';

function Selected(props) {
  return (
    <Row>
      <Col sm={12}>
        <h1>Selected in Map</h1>
        <hr/>
        <MolsByMolsets
          {...props}
          component={MolsByMolsetsTabs}
        />
      </Col>
    </Row>
  )
}

function MapTabs(props) {
  const [selectedMolsInMap, setSelectedMolsInMap] = React.useState([]);
  const [selectedPointsInMap, setSelectedPointsInMap] = React.useState([]);
  const [activeTab, setActiveTab] = React.useState('Map');

  const tabs = {
    map: {
      title: "Map"
    },
    selected: {
      title: "Selected"
    }
  };

  return (
    // <TabWidgetSmart
    //   {...props}
    //   tabs={tabs}
    //   selectedMols={selectedMolsInMap}
    //   selectedPoints={selectedPointsInMap}
    //   onMolsSelect={setSelectedMolsInMap}
    //   onPointsSelect={setSelectedPointsInMap}
    // />
    <Tabs
      defaultActiveKey={tabs.map.title}
      activeKey={activeTab}
      onSelect={k => setActiveTab(k)}
      id="noanim-tab-example"
    >
      <Tab
        eventKey={tabs.map.title}
        title={tabs.map.title}
        // className={classnames({ active: activeTab === tabs.map.title })}
      >
        <MapPage
          {...props}
          selectedMols={selectedMolsInMap}
          selectedPoints={selectedPointsInMap}
          onMolsSelect={setSelectedMolsInMap}
          onPointsSelect={setSelectedPointsInMap}
        />
      </Tab>

      <Tab
        eventKey={tabs.selected.title}
        title={tabs.selected.title}
        // className={classnames({ active: activeTab === tabs.selected.title })}
      >
        <Selected
          {...props}
          selectedMols={selectedMolsInMap}
          selectedPoints={selectedPointsInMap}
          onMolsSelect={setSelectedMolsInMap}
          onPointsSelect={setSelectedPointsInMap}
        />
      </Tab>

    </Tabs>
  )
}

const MapDashboard = (props) => {
  const [selectedMap, setSelectedMap] = React.useState(null);

  const defaultMapClass = "Map";
  return (
    <ComponentWithObjects
      objectListURL={props.apiUrls.mapsRoot}
      emptyClassName={defaultMapClass}
      currentProject={props.currentProject}
      render={
        (mapObjects) => {
          const maps = mapObjects[defaultMapClass];
          const selected = selectedMap ? selectedMap : (
            maps.length > 0 ? setSelectedMap(maps[0]) : null
          );
          if (selected) {
            const molsets = selected.molsets;
            const resourcesDef = {};
            molsets.forEach(molset => {
              molset.activities.forEach(activitySetID => {
                resourcesDef[activitySetID] = new URL(`${activitySetID}/`, props.apiUrls.activitySetsRoot)
              })
            });
            return (
              <ComponentWithResources
                definition={resourcesDef}
              >
                {
                  (allLoaded, activitySets) => {
                    return allLoaded ? (
                      <MapTabs
                        {...props}
                        selectedMap={selected}
                        maps={maps}
                        molsets={molsets}
                        activitySets={activitySets}
                        onMapSelect={setSelectedMap}
                      />
                    ) : <div>Loading...</div>
                  }
                }
              </ComponentWithResources>
            )
          } else {
            return <div>Select a map from the menu above.</div>
          }
        }
      }
    />
  )
};

export default MapDashboard;
