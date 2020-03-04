import React from "react";
import { CardSubtitle } from 'reactstrap';

function ProjectItemSubTitle(props) {
  const item = props.item;
  return (
    <React.Fragment>
      <CardSubtitle>
        <p>
          Created: {
          new Date(item.created).toLocaleDateString()
          + ' – ' + new Date(item.created).toLocaleTimeString()
        }
          <br/>
          Last Update: {
          new Date(item.updated).toLocaleDateString()
          + ' – ' + new Date(item.updated).toLocaleTimeString()
        }
        </p>
      </CardSubtitle>
    </React.Fragment>
  )
}

export default ProjectItemSubTitle;