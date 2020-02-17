import React from "react";
import { ComponentWithResources, ResponsiveGrid, TaskAwareComponent } from '../../index';
import { Card } from 'reactstrap';

class ModelGrid extends React.Component {

  render() {
    const chosenAlgorithm = this.props.chosenAlgorithm;
    const models = this.props.models;
    const newCardSetup = this.props.newCardSetup ? this.props.newCardSetup : {
      h : {"md" : 15, "sm" : 15},
      w : {"md" : 1, "sm" : 1},
      minH : {"md" : 3, "sm" : 3},
    };
    const cardSetup = this.props.cardSetup ? this.props.cardSetup : {
      h : {"md" : 12, "sm" : 12},
      w : {"md" : 1, "sm" : 1},
      minH : {"md" : 3, "sm" : 3},
    };

    if (models.length === 0 && !chosenAlgorithm) {
      return <p>Start by selecting an algorithm. See the actions menu in the top right.</p>
    }

    const existing_cards = models.map(model => (Object.assign({
      id : model.id,
      data : model
    }, cardSetup)));
    const new_card = Object.assign({
      id : "new-model",
    }, newCardSetup);

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
                        // <LiveObject {...this.props} url={new URL(`${item.data.id}/`, this.props.listURL)}>
                        <ComponentWithResources
                          definition={{
                            model : new URL(`${item.data.id}/`, this.props.listURL),
                          }}
                          updateAfterTasksDone={true}
                          {...taskInfo}
                        >
                        {
                            (loaded, resources) => loaded ? (
                              <ModelComponent
                                {...this.props}
                                {...taskInfo}
                                onTaskUpdate={onTaskUpdate}
                                model={resources.model}
                                modelUrl={new URL(`${item.data.id}/`, this.props.listURL)}
                                modelClass={this.props.modelClass}
                                onModelDelete={this.props.handleModelDelete}
                              />
                            ) : <div>Loading...</div>
                          }
                        </ComponentWithResources>
                        // </LiveObject>
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