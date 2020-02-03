import React from "react";
import { LiveObject, ResponsiveGrid, TaskAwareComponent } from '../../index';
import { Card } from 'reactstrap';

class ModelGrid extends React.Component {
  render() {
    const chosenAlgorithm = this.props.chosenAlgorithm;
    if (Object.entries(this.props.models).length === 0 && !chosenAlgorithm) {
      return <p>Start by selecting a QSAR modelling algorithm in the actions menu.</p>
    }

    const models = this.props.models[this.props.modelClass] ? this.props.models[this.props.modelClass] : [];
    const existing_cards = models.map(model => ({
      id : model.id,
      h : {"md" : 12, "sm" : 12},
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

    const ModelComponent = this.props.modelComponent;
    const NewModelComponent = this.props.newModelComponent;

    return (
      <div className="models-grid">
        <ResponsiveGrid
          items={existing_cards.concat(new_card)}
          rowHeight={75}
          mdCols={2}
          smCols={1}
        >
          {
            existing_cards.map(
              item => (
                <Card key={item.id}>
                  <TaskAwareComponent
                    handleResponseErrors={this.props.handleResponseErrors}
                    tasksURL={new URL(`${item.data.id}/tasks/all/`, this.props.listURL)}
                    render={
                      (taskInfo, onTaskUpdate) => (
                        <LiveObject {...this.props} url={new URL(`${item.data.id}/`, this.props.listURL)}>
                          {
                            (model) => (
                              <ModelComponent
                                {...this.props}
                                {...taskInfo}
                                onTaskUpdate={onTaskUpdate}
                                model={model}
                                modelClass={this.props.modelClass}
                                onModelDelete={this.props.handleModelDelete}
                              />
                            )
                          }
                        </LiveObject>
                      )
                    }
                  />
                </Card>
              )
            ).concat(chosenAlgorithm ? [(
              <Card key={new_card.id} id={new_card.id}>
                <NewModelComponent {...this.props}/>
              </Card>
            )] : [])
          }
        </ResponsiveGrid>
      </div>
    )
  }
}

export default ModelGrid