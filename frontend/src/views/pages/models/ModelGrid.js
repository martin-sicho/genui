import React from "react";
import { ResponsiveGrid, TaskAwareComponent } from '../../../genui';
import { Card } from 'reactstrap';
import ModelCard from './ModelCard';
import ModelCardNew from './ModelCardNew';

class ModelGrid extends React.Component {
  MODEL_CLASS = "QSARModel";

  render() {
    const chosenAlgorithm = this.props.chosenAlgorithm;
    if (Object.entries(this.props.models).length === 0 && !chosenAlgorithm) {
      return <p>Start by selecting a QSAR modelling algorithm in the actions menu.</p>
    }

    const models = this.props.models[this.MODEL_CLASS] ? this.props.models[this.MODEL_CLASS] : [];
    const existing_cards = models.map(model => ({
      id : model.id,
      h : {"md" : 9, "sm" : 8},
      w : {"md" : 1, "sm" : 1},
      minH : {"md" : 3, "sm" : 3},
      data : model
    }));
    const new_card = {
      id : "new-model",
      h : {"md" : 15, "sm" : 15},
      w : {"md" : 1, "sm" : 1},
      minH : {"md" : 3, "sm" : 3},
    };

    return (
      <ResponsiveGrid
          items={existing_cards.concat(new_card)}
          rowHeight={75}
          mdCols={2}
          smCols={1}
        >
          {
            existing_cards.map(
              item => (
                <Card key={item.id.toString()}>
                  <TaskAwareComponent
                    handleResponseErrors={this.props.handleResponseErrors}
                    tasksURL={new URL(`models/${item.data.id}/tasks/all/`, this.props.apiUrls.qsarRoot)}
                    render={
                      (taskInfo, onTaskUpdate) => (
                        <ModelCard
                          {...this.props}
                          {...taskInfo}
                          onTaskUpdate={onTaskUpdate}
                          model={item.data}
                          modelClass={this.MODEL_CLASS}
                          onModelDelete={this.props.handleModelDelete}
                        />
                      )
                    }
                  />
                </Card>
              )
            ).concat(chosenAlgorithm ? [(
              <Card key={new_card.id} id={new_card.id}>
                <ModelCardNew {...this.props}/>
              </Card>
            )] : [])
          }
        </ResponsiveGrid>
    )
  }
}

export default ModelGrid