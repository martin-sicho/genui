import React from 'react';
import { Col, Row } from 'reactstrap';
import {
  ActivitiesByTypeFlatView,
  MoleculeActivityProvider,
  MoleculeMetadata,
  MoleculeImage,
  MoleculePropsProvider,
  TabWidget, PropertiesTable,
} from '../../..';
import SimplePaginator from '../../SimplePaginator';

class MoleculeData extends React.Component {
  constructor(props) {
    super(props);

    this.state = this.initState(props);
  }

  initState = (props) => {
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
            component={ActivitiesByTypeFlatView}
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

    return {
      showData : showData,
      showActivities: showActivities,
      showProperties: showProperties,
      tabs: tabs
    };
  };

  shouldComponentUpdate(nextProps, nextState, nextContext) {
    if (this.props.updateCondition) {
      return this.props.updateCondition(this.props, nextProps, this.state, nextState, nextContext);
    } else {
      return true;
    }
  }

  render() {
    return (
      <TabWidget {...this.props} tabs={this.state.tabs} activeTab={this.state.showActivities ? "Activities" : "Info"}/>
    )
  }
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
        <MoleculeData {...props}/>
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