import {
  Card,
  CardBody,
  Col,
  Row,
} from 'reactstrap';
import React from 'react';
import Map from './Map';
import { MolsByMolsets, MolsByMolsetsTabs } from './MolsByMolsets';
import { MoleculeDetail, MoleculeActivityDetail } from '../../../../genui';
import MapHandlers from './MapHandlers';
import MapSelectorPage from './MapSelectorPage';

function ActivitiesList(props) {
  const set = props.set;
  const activities = props.activities;

  return (
    <React.Fragment>
      <h6>{set.name}</h6>
      <ul>
        {
          activities.map(item => {
            return (
              <li key={item.id}>{item.value}</li>
            )
          })
        }
      </ul>
    </React.Fragment>
  )
}

function ActivitySetList(props) {
  const actSets = props.activitySets;
  const activities = props.activities;

  return (
    <ul>
      {
        Object.keys(actSets).map(key => {
          const set = actSets[key];
          if (activities[key].length > 0) {
            return <li key={set.id}><ActivitiesList set={set} activities={activities[key]}/></li>
          } else {
            return null;
          }
        })
      }
    </ul>
  )
}

class MapsPageComponents extends React.Component {

  render() {
    const hoverMol = this.props.hoverMol;
    const selectedMap = this.props.map;
    return (
      selectedMap ? (

        <React.Fragment>
          <Row>
            <Col md={8} sm={10}>
              <Card>
                <CardBody>
                  <Map
                    {...this.props}
                  />
                </CardBody>
              </Card>

            </Col>

            <Col md={4} sm={2}>
              {
                hoverMol ? (
                  <Row>
                    <Col md={6} sm={4}>
                      <MoleculeDetail
                        {...this.props}
                        mol={hoverMol}
                      />
                    </Col>
                    <Col md={6} sm={8}>
                      <MoleculeActivityDetail
                        {...this.props}
                        mol={hoverMol}
                        component={ActivitySetList}
                      />
                    </Col>
                  </Row>
                ) : null
              }

              <hr/>

              <div>Activity and summary stats for the displayed molecule sets...</div>
            </Col>
          </Row>
          <Row>
            <Col sm={12}>
              <h1>Selected Compounds</h1>
              <hr/>
              <MolsByMolsets
                {...this.props}
                component={MolsByMolsetsTabs}
              />
            </Col>
          </Row>
        </React.Fragment>

      ) : <div>Select a map to display from the menu.</div>
    )
  }
}

function MapsPage(props) {

  const PageImpl = props => {
    return <MapHandlers {...props} component={MapsPageComponents}/>
  };
  return (
    <MapSelectorPage {...props} component={PageImpl}/>
    )
}

export default MapsPage;