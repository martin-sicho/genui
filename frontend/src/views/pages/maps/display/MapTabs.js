import React from 'react';
import MapPage from './MapPage';
import { TabWidget} from '../../../../genui';
import SelectedActivitiesPage from './SelectedActivitiesPage';
import SelectedList from './SelectedList';
import { Col, Row } from 'reactstrap';

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
      title: "List",
      renderedComponent: (props) => (
        <React.Fragment>
          <h1>List of Selected Compounds in the Map</h1>
          <hr/>
          <Row>
            <Col sm={12}>
              {
                props.selectedMols.length > 0 ? <SelectedList {...props}/> : <p>Select compounds in the map to see details.</p>
              }
            </Col>
          </Row>
        </React.Fragment>
      )
    },
    {
      title: "Activities",
      renderedComponent: SelectedActivitiesPage
    },
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