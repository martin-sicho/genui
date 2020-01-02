import React from 'react';
import { Col, Row } from 'reactstrap';

class MolsStats extends React.Component {

  constructor(props) {
    super(props);

    this.state = {
      molCount : 0
    }
  }

  componentDidMount() {
    fetch(this.props.moleculesURL)
      .then(response => response.json())
      .then(this.updateMolStats)
  }

  updateMolStats = (data) => {
    this.setState({molCount : data.count});
  };

  render() {
    return (
      <React.Fragment>
        <h4>Compounds</h4>
        <p>Total: {this.state.molCount}</p>
        <h4>Associated Targets</h4>
        <ul>
          {
            this.props.molset.targets.map(
              target => <li key={target.targetID}>{target.targetID}</li>
            )
          }
        </ul>
      </React.Fragment>
    )
  }
}

class MolSetImportStatus extends React.Component {

  render() {
    return (
      <React.Fragment>
        <h4>Started Tasks</h4>
        <p>Something...</p>
        <h4>Completed Tasks</h4>
        <p>Something...</p>
      </React.Fragment>
    )
  }
}

class ChEMBLInfo extends React.Component {

  constructor(props) {
    super(props);

    this.molset = this.props.molset;
  }

  render() {
    return (
      <Row>
        <Col sm="12">
          <h4>Description</h4>
          <p>{this.molset.description}</p>
          <MolsStats molset={this.props.molset} moleculesURL={this.props.moleculesURL} />
          <MolSetImportStatus molset={this.props.molset} tasksURL={this.props.tasksURL}/>
        </Col>
      </Row>
    );
  }
}

export default ChEMBLInfo;