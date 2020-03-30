import React from "react"
import { Button, Col, Row, Table } from 'reactstrap';
import { ComponentWithObjects } from '../../../../index';
import PredictionsCreateForm from './PredictionsCreateForm';

function PredsList(props) {
  const activitySets = props.activitySets;
  const molsets = props.molsets;
  return (
    activitySets.length > 0 ? (
      <div className="model-predictions-list">
      <Table>
        <thead>
          <tr>
            <th>Name</th>
            <th>Compounds</th>
            <th>Delete</th>
          </tr>
        </thead>

        <tbody>
          {
            activitySets.map(set => {
              const molset = molsets.find(item => item.id === set.molecules);
              return (
                <tr key={set.id}>
                  <th scope="row">{set.name}</th>
                  <td>{molset.name}</td>
                  <td>
                    <Button color="danger" onClick={() => props.handleDeleteActivitySet(set.className, set)}> <i className="fa fa-close"/>Delete</Button>
                  </td>
                </tr>
              )
            })
          }
        </tbody>
      </Table>

    </div>
    ) : <p>No predictions, yet. You can create one below.</p>
  )
}

function ModelPredsPage(props) {
  return (
    <React.Fragment>
      <Row>
        <Col sm="12">
          <h4>Existing Predictions</h4>
          <PredsList {...props}/>
        </Col>
      </Row>
      <Row>
        <Col sm="12">
          <h4>Create New Predictions</h4>
          <PredictionsCreateForm {...props}/>
        </Col>
      </Row>
    </React.Fragment>
  )
}

class ModelPreds extends React.Component {

  render() {
    let molsets = [];
    Object.keys(this.props.compoundSets).forEach(csetClass => {
      molsets = molsets.concat(this.props.compoundSets[csetClass]);
    });
    const defaultClass = "ModelActivitySet";
    const listUrl = new URL(`models/${this.props.model.id}/predictions/`, this.props.apiUrls.qsarRoot);
    return (
      <ComponentWithObjects
        objectListURL={listUrl}
        emptyClassName={defaultClass}
        currentProject={this.props.currentProject}
        customDelete={(className, toDelete) => {
          const url = new URL(`${toDelete.id}/`, this.props.apiUrls.activitySetsRoot);
          fetch(url, {method: 'DELETE', credentials: "include",})
            .catch(
              (error) => console.log(error)
            );
        }}
        render={
          (
            activitySets,
            handleAddActivitySetList,
            handleAddActivitySet,
            handleDeleteActivitySet,
          ) => {
            return (<ModelPredsPage
              {...this.props}
              molsets={molsets}
              predictionsListUrl={listUrl}
              defaultClass={defaultClass}
              activitySets={activitySets[defaultClass]}
              handleAddActivitySet={handleAddActivitySet}
              handleDeleteActivitySet={handleDeleteActivitySet}
            />)
          }
        }
      />
    )
  }
}

export default ModelPreds;