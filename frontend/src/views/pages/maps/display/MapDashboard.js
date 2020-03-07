import React from 'react';
import { ComponentWithObjects, ComponentWithResources } from '../../../../genui';
import MapPage from './MapPage';
import { Col, Row } from 'reactstrap';
import { MolsByMolsets, MolsByMolsetsTabs } from './MolsByMolsets';
import { Tab, Tabs } from 'react-bootstrap';
import HeaderNav from './HeaderNav';

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

class MapDashboard extends React.Component {

  constructor(props) {
    super(props);

    this.state ={
      selectedMap: null
    }
  }

  componentDidUpdate(prevProps, prevState, snapshot) {
    const maps = this.props.maps;
    const selectedMap = this.state.selectedMap;
    if (maps.length > 0 && !selectedMap) {
      this.props.onHeaderChange(
        <HeaderNav
          {...this.props}
          maps={maps}
          onMapChoice={map => {
            this.setState({
              selectedMap: map,
            })
          }}/>
      );

      this.setState({
        selectedMap: maps[0]
      })
    }
  }

  render() {
    const maps = this.props.maps;
    const selected = this.state.selectedMap;

    if (selected) {
      const molsets = selected.molsets;
      const resourcesDef = {};
      molsets.forEach(molset => {
        molset.activities.forEach(activitySetID => {
          resourcesDef[activitySetID] = new URL(`${activitySetID}/`, this.props.apiUrls.activitySetsRoot)
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
                  {...this.props}
                  selectedMap={selected}
                  maps={maps}
                  molsets={molsets}
                  activitySets={activitySets}
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

function Maps(props) {
  const defaultMapClass = "Map";
  return (
    <ComponentWithObjects
      objectListURL={props.apiUrls.mapsRoot}
      emptyClassName={defaultMapClass}
      currentProject={props.currentProject}
      render={
        (mapObjects) => {
          const maps = mapObjects[defaultMapClass];
          return <MapDashboard {...props} maps={maps}/>
        }
      }
    />
  )
}

export default Maps;
