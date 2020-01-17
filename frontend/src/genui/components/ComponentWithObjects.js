import React from "react";

class ComponentWithObjects extends React.Component {

  constructor(props) {
    super(props);

    this.objectListRoot = this.props.objectListURL;
    this.emptyClassProperty = props.classProperty ? props.classProperty : "no_class_property";

    this.state = {
      fetchUpdates : false,
      objects : null,
    }
  }

  componentDidMount() {
    this.setState({fetchUpdates : true})
  }

  componentDidUpdate(prevProps, prevState, snapshot) {
    if (this.props.currentProject && this.state.fetchUpdates) {
      this.fetchUpdates(this.props.currentProject);
    }
  }

  fetchUpdates = (project) => {
    const params = new URLSearchParams();
    params.append('project_id', project.id);
    fetch(this.objectListRoot.toString() + "?" + params.toString())
      .then(response => response.json())
      .then(this.getObjects)
  };

  getObjects = (data) => {
    const objects = {};
    for (const obj of data) {
      if (!obj.className) {
        obj.className = this.emptyClassProperty
      }
      if (!objects.hasOwnProperty(obj.className)) {
        objects[obj.className] = [];
      }
      objects[obj.className].push(obj);
    }
    this.setState({
      objects : objects,
      fetchUpdates : false
    })
  };

  handleAddObject = (className, data) => {
    this.setState(prevState => {
      prevState.objects[className].push(data);
      return {
        objects : prevState.objects
      };
    });
  };

  handleAddObjectList = (className, objectList, overwrite=false) => {
    this.setState((prevState) => {
      const old_objects = prevState.objects;
      if (old_objects.hasOwnProperty(className)) {
        if (overwrite) {
          old_objects[className] = objectList;
        } else {
          old_objects[className] = old_objects[className].concat(objectList);
        }
      } else {
        old_objects[className] = objectList;
      }
      return {
        objects : old_objects
      }
    });
  };

  handleObjectDelete = (className, object) => {
    fetch(this.objectListRoot.toString() + object.id + '/', {method: 'DELETE'})
      .then(
        () => {
          this.setState(prevState => {
            const object = prevState.objects[className];
            const idx_del = object.findIndex(item => item.id === object.id);
            object.splice(idx_del, 1);
            return {
              objects : prevState.objects
            };
          });
        }
      ).catch(
      (error) => console.log(error)
    )
    ;
  };

  render() {
    if (this.state.objects === null) {
      return <div>Loading...</div>
    }

    return (
      <React.Fragment>
        {this.props.render(this.state.objects, this.handleAddObjectList, this.handleAddObject, this.handleObjectDelete)}
      </React.Fragment>
    )
  }
}

export default ComponentWithObjects;