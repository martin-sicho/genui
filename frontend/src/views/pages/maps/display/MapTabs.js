import React from 'react';
import { Tab, Tabs } from 'react-bootstrap';
import MapPage from './MapPage';
import { Col, Row } from 'reactstrap';
import {MolsToMolSetGroups, MolSetsTabs} from '../../../../genui';

function Selected(props) {
  return (
    <Row>
      <Col sm={12}>
        <h1>Selected in Map</h1>
        <hr/>
        <MolsToMolSetGroups
          {...props}
          mols={props.selectedMols}
        >
          {
            (groups) => {
              return (
                <MolSetsTabs
                  {...props}
                  groupedMols={groups}
                  mols={props.selectedMols}
                />
              )
            }
          }
        </MolsToMolSetGroups>
      </Col>
    </Row>
  )
}

export default function MapTabs(props) {
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
    <Tabs
      defaultActiveKey={tabs.map.title}
      activeKey={activeTab}
      onSelect={k => setActiveTab(k)}
      id="noanim-tab-example"
    >
      <Tab
        eventKey={tabs.map.title}
        title={tabs.map.title}
        // className={classnames({ active: activeTab === summaries.map.title })}
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
        // className={classnames({ active: activeTab === summaries.selected.title })}
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