import React from 'react';
import { Card, CardTitle, Col, Row } from 'reactstrap';
import { MoleculeImage, MoleculePropsProvider, PropertiesTable } from '../../../../genui';
import ActivitySummaryPlotter from './ActivitySummaryPlotter';
import SelectedListPage from './SelectedListPage';

function ActivitySummaryPlot(props) {
  const [hoverMol, setHoverMol] = React.useState(null);

  return (
    <Row>
      <Col sm={8}>
        <ActivitySummaryPlotter
          {...props}
          mols={props.selectedMols}
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

export default function ActivitySummary(props) {

  const [selectedInOverview, setSelectedInOverview] = React.useState([]);
  const [selectedInOverviewRev, setSelectedInOverviewRev] = React.useState(0);

  return (
    <React.Fragment>
      <ActivitySummaryPlot
        {...props}
        onMolsSelect={setSelectedInOverview}
        setSelectedMolsOverviewRevision={setSelectedInOverviewRev}
        selectedMolsOverviewRevision={selectedInOverview}
      />
      <hr/>
      <SelectedListPage
        {...props}
        selectedMols={selectedInOverview}
        selectedMolsRevision={selectedInOverviewRev}
        emptyMessage="Select points in the plot above to see the list of compounds associated with those activities."
      />
    </React.Fragment>
  )
}