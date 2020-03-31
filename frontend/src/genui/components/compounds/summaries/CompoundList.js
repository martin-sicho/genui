import React from 'react';
import { Col, Row } from 'reactstrap';
import {
  ActivitySetTabView,
  MoleculeActivityProvider,
  MoleculeMetadata,
  MoleculeImage,
  MoleculePropsProvider,
  TabWidget, PropertiesTable,
} from '../../..';
import SimplePaginator from '../../SimplePaginator';

function DataTabs(props) {
  let showData = typeof(props.showInfo) === 'boolean' ? props.showInfo : true;
  const showActivities = typeof(props.showActivities) === 'boolean' ? props.showActivities : true;
  const showProperties = typeof(props.showProperties) === 'boolean' ? props.showProperties : true;

  if (!(showData || showActivities)) {
    showData = true;
  }

  const tabs = [];
  if (showData) {
    tabs.push({
      title: "Info",
      renderedComponent: MoleculeMetadata
    },)
  }

  if (showActivities) {
    tabs.push({
      title: "Activities",
      renderedComponent: (props) => (
        <MoleculeActivityProvider
          {...props}
          component={ActivitySetTabView}
        />
      )
    })
  }

  if (showProperties) {
    tabs.push({
      title: "Properties",
      renderedComponent: (props) => (
        <MoleculePropsProvider
          {...props}
          propsList={[
            "AMW",
            "NUMHEAVYATOMS",
            "NUMAROMATICRINGS",
            "HBA",
            "HBD",
            "LOGP",
            "TPSA",
          ]}
          component={PropertiesTable}
        />
      )
    })
  }

  return (
    <TabWidget {...props} tabs={tabs} activeTab={showData ? "Info" : "Activities"}/>
  )
}

export function CompoundListItem(props) {
  const mol = props.mol;
  const sm_cols = [3, 9];
  const md_cols = [3, 9];

  return (
    <Row>
      <Col md={md_cols[0]} sm={sm_cols[0]}>
        <MoleculeImage mol={mol}/>
      </Col>
      <Col md={md_cols[1]} sm={sm_cols[1]}>
        <DataTabs {...props}/>
      </Col>
    </Row>
  )
}

export function CompoundListPageItem(props) {
  return <CompoundListItem {...props} mol={props.pageItem}/>
}

export default function CompoundList(props) {
  const mols = props.mols;

  if (props.paginate) {
    return (
      <SimplePaginator {...props} items={mols} component={CompoundListPageItem}/>
    )
  } else {
      return (
        <React.Fragment>
          {
            mols.map(mol => (
              <CompoundListItem {...props} key={mol.id} mol={mol}/>
            ))
          }
        </React.Fragment>
      )
  }
}