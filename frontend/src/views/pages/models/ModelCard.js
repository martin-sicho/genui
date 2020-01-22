import React from "react";
import { Button, CardBody, CardFooter, CardHeader } from 'reactstrap';
import { TabWidget, ProjectItemSubTitle } from '../../../genui';
import ModelInfo from './tabs/ModelInfo';
import ModelPerformance from './tabs/ModelPerformance';
import ModelPredictions from './tabs/ModelPredictions';

class ModelCard extends React.Component {

  constructor(props) {
    super(props);

    this.state = {
      isDeleting : false
    }
  }

  render() {
    const model = this.props.model;
    const tabs = [
      {
        title : "Info",
        renderedComponent : () =>
          <ModelInfo
            {...this.props}
          />
      },
      {
        title: "Performance"
        , renderedComponent : () =>
          <ModelPerformance
            {...this.props}
          />
      },
      {
        title: "Predictions"
        , renderedComponent : () =>
          <ModelPredictions
            {...this.props}
          />
      }
    ];

    return (
      <React.Fragment>
        <CardHeader>{model.name}</CardHeader>

        <CardBody className="scrollable">
          <ProjectItemSubTitle item={model}/>
          <TabWidget tabs={tabs} />
        </CardBody>

        <CardFooter>
          <Button color="danger" disabled={this.state.isDeleting} onClick={() => {this.setState({isDeleting:true});this.props.onModelDelete(this.props.modelClass, model)}}>Delete</Button>
        </CardFooter>
      </React.Fragment>
    )
  }
}

export default ModelCard;