import React from 'react';
import { Col, Row } from 'reactstrap';
import Map from './Map';
import MapSidebar from './MapSidebar';

const MapDashboard = (props) => {
  return (
    <Row>
      <Col md={8} sm={10}>
        <Map {...props}/>
      </Col>

      <Col md={4} sm={2}>
        <MapSidebar {...props}/>
      </Col>
    </Row>
  );
};

export default MapDashboard;
