import React from 'react';
import { Card, CardTitle, Col, Row } from 'reactstrap';
import { MoleculeImage, MoleculePropsProvider, PropertiesTable } from '../../../../genui';
import ActivitySummaryPlotter from './ActivitySummaryPlotter';
import SelectedList from './SelectedList';

function PlotComponent(props) {
  const [hoverMol, setHoverMol] = React.useState(null);

  return (
    <Row>
      <Col sm={8}>
        <ActivitySummaryPlotter
          {...props}
          onMolHover={setHoverMol}
        />
      </Col>
      <Col sm={4}>
        {hoverMol ? (
          <React.Fragment>
            <MoleculeImage mol={hoverMol}/>
            <hr/>
            <Card body>
              <CardTitle>Properties</CardTitle>
              {/*<CardText>*/}
              {/*<MoleculeMetadata mol={hoverMol}/>*/}
              {/*</CardText>*/}
              <MoleculePropsProvider
                {...props}
                mol={hoverMol}
                propsList={[
                  "AMW",
                  "NUMHEAVYATOMS",
                  "NUMAROMATICRINGS",
                  "HBA",
                  "HBD",
                  "LOGP",
                  "TPSA",
                ]}
                updateCondition={(prevProps, currentProps) => {
                  return prevProps.mol && (prevProps.mol.id !== currentProps.mol.id)
                }}
                component={PropertiesTable}
              />
            </Card>
          </React.Fragment>
        ) : <p>Hover over a point in the plot to see details.</p>}
      </Col>
    </Row>
  )
}

function ListComponent(props) {
  const selectedInOverview = props.selectedInOverview;
  return (
    <Row>
      <Col sm={12}>
        <h2>Selected Compounds</h2>
        {
          selectedInOverview.length > 0 ? (
            <SelectedList
              {...props}
              selectedMols={selectedInOverview}
            />
          ) : <p>Select points in the summary plot above to see details.</p>
        }
      </Col>
    </Row>
  )
}

export default function ActivitySummary(props) {

  const [selectedInOverview, setSelectedInOverview] = React.useState([]);

  return (
    <React.Fragment>
      <PlotComponent {...props} onMolsSelect={setSelectedInOverview}/>
      <hr/>
      <ListComponent {...props} selectedInOverview={selectedInOverview}/>
    </React.Fragment>
  )
}